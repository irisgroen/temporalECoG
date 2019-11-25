function [err, pred] = DN(param, data, stim, srate)
%
% function [err, pred] = DN(param, data, stim, srate)
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
x       = []; % a struct of model parameteres

numtimepts  = size(stim,1);
numstim     = size(stim,2);

normSum = @(x) x./sum(x);

%% SET UP THE MODEL PARAMETERS

fields = {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};
x      = toSetField(x, fields, param);

%% COMPUTE THE IMPULSE RESPONSE FUNCTION

% HERE I ASSUME THAT THE NEGATIVE PART OF THE IMPULSE RESPONSE HAS A TIME
% CONSTANT 1.5 TIMES THAT OF THE POSITIVE PART OF THE IMPULSE RESPONSE
%if x.tau1 > 0.5, warning('tau1>1, the estimation for other parameters may not be accurate'); end
    
t       = (1:numtimepts)' / srate;

irf_pos = gammaPDF(t, x.tau1, 2);
irf_neg = gammaPDF(t, x.tau1*1.5, 2);
irf     = irf_pos - x.weight.* irf_neg;

%% COMPUTE THE DELAYED REPSONSE FOR THE NORMALIZATION POOL

irf_norm = normSum(exp(-t/x.tau2));

%% COMPUTE THE NORMALIZATION RESPONSE

% ADD SHIFT TO THE STIMULUS -------------------------------------------
sft       = round(x.shift * srate);
stimtmp   = padarray(stim, [sft, 0], 0, 'pre');
stim = stimtmp(1 : size(stim, 1), :);

% COMPUTE THE NORMALIZATION NUMERATOR ---------------------------------
linrsp  = conv2(stim, irf, 'full'); linrsp = linrsp(1:numtimepts,:);
numrsp  = abs(linrsp).^x.n;

% COMPUTE THE NORMALIZATION DENOMINATOR -------------------------------
poolrsp = conv2(linrsp, irf_norm, 'full'); poolrsp = poolrsp(1:numtimepts,:);
demrsp  = x.sigma.^x.n + abs(poolrsp).^x.n;

% COMPUTE THE NORMALIZATION RESPONSE
normrsp = x.scale.*(numrsp./demrsp);


pred = normrsp;

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end
end