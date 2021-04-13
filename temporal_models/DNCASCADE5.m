function [err, pred] = DNCASCADE(param, data, stim, srate)
%
% function [err, pred] = DNCASCADE(param, data, stim, srate)
% INPUTS  -----------------------------------------------------------------
% params : 7 fields.
%          1. tau1 -- time to peak for positive IRF (seconds)
%          2. weight -- ratio of negative to positive IRFs (set to 0
%               for uniphasic irf)
%          3. tau2 -- time window of adaptation (seconds)
%          4. n -- exponent
%          5. sigma -- semi-saturation constant
%          6. shift -- time between stimulus onset and when the signal reaches
%               the cortex (seconds)
%          7. scale -- response gain.
%
% data :   matrix, samples x trials
%
% stim :   matrix, samples x trial 
%
% srate :  sample rate in Hz
%
% OUTPUTS -----------------------------------------------------------------
% err:  sum of squared error
% pred: predicted time series


%% PRE-DEFINED /EXTRACTED VARIABLES

numtimepts  = size(stim,1);

normSum = @(x) x./sum(x);

%% SET UP THE MODEL PARAMETERS

fields = {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};
prm      = toSetField([], fields, param);

%% COMPUTE THE IMPULSE RESPONSE FUNCTION   
t       = (1:numtimepts)' / srate;

irf_pos = gammaPDF(t, prm.tau1, 2);
irf_neg = gammaPDF(t, prm.tau1*1.5, 2);
irf     = irf_pos - prm.weight.* irf_neg;

%% COMPUTE THE DELAYED REPSONSE FOR THE NORMALIZATION POOL

irf_norm = normSum(exp(-t/prm.tau2));

%% COMPUTE THE NORMALIZATION RESPONSE


% ADD SHIFT TO THE STIMULUS -------------------------------------------
t2 = t-prm.shift;
stim = interp1(t, stim, t2, [], 0);

ncascades = 5;

rsp = stim;
for ii = 1:ncascades
    rsp = dncomputeOneLayer(rsp, irf, irf_norm, prm);
end

pred = rsp;

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end

end

function rsp = dncomputeOneLayer(stim, irf, irf_norm, x)

n = size(stim,1);

% COMPUTE THE NORMALIZATION NUMERATOR ---------------------------------
linrsp  = conv2(stim, irf, 'full'); linrsp = linrsp(1:n,:);
numrsp  = abs(linrsp).^x.n;

% COMPUTE THE NORMALIZATION DENOMINATOR -------------------------------
poolrsp = conv2(linrsp, irf_norm, 'full'); poolrsp = poolrsp(1:n,:);
demrsp  = x.sigma.^x.n + abs(poolrsp).^x.n;

% COMPUTE THE NORMALIZATION RESPONSE
rsp = x.scale.*(numrsp./demrsp);

end