function [err, pred] = LINEAR_RECTF_EXP_NORM_DELAY(param, data, stim, srate)
%
% function [err, pred, numrsp, demrsp] = LINEAR_RECTH_EXP_NORM_DELAY(param, data, stim, srate)
% INPUTS  -----------------------------------------------------------------
% params : 9 fields.
%          1. tau1 -- time to peak for positive IRF (seconds)
%          2. tau2 -- time to peak for negative IRF (seconds)
%          3. n_irf -- phase delay for positive and negative IRF 
%          4. weight -- ratio of negative to positive IRFs 
%          5. shift -- time between stimulus onset and when the signal reaches
%               the cortex (seconds)
%          6. scale -- response gain.
%          7. n -- exponent 
%          8. sigma -- semi-saturation constant
%          9. tau_a -- time window of adaptation (seconds)
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

%% SET UP THE MODEL PARAMETERS
fields    = {'tau1', 'tau2', 'n_irf', 'weight', 'shift', 'scale', 'n', 'sigma', 'tau_a'};
prm       = toSetField([], fields, param);
prm.n_irf = max(round(prm.n_irf),1); % n_irf has to be an integer and can't be zero

%% COMPUTE THE IMPULSE RESPONSE FUNCTION

t       = (1:numtimepts)' / srate;

% COMPUTE THE IRF
irf_pos = gammaPDF(t, prm.tau1, prm.n_irf);
irf_neg = gammaPDF(t, prm.tau2, prm.n_irf);
irf     = irf_pos - prm.weight.* irf_neg;

%% CREATE EXPONENTIAL DECAY FILTER

normSum = @(x) x./sum(x);
irf_norm = normSum(exp(-t/prm.tau_a));

%% COMPUTE THE RESPONSE

% ADD SHIFT TO THE STIMULUS -------------------------------------------
t2 = t-prm.shift;
stim = interp1(t, stim, t2, [], 0);

% COMPUTE THE CONVOLUTION
linrsp  = conv2(stim, irf, 'full');         % convolve
linrsp  = linrsp(1:numtimepts,:);           % cut
numrsp  = abs(linrsp);                      % full wave rectification
numrsp  = numrsp.^prm.n;                    % exponentiate

% COMPUTE THE NORMALIZED RESPONSE
poolrsp = conv2(linrsp, irf_norm, 'full');        % convolve
poolrsp = poolrsp(1:numtimepts,:);                % cut
demrsp  = prm.sigma.^prm.n + abs(poolrsp).^prm.n; % semi-sat + exponentiate
normrsp = numrsp./demrsp;                         % divide 

% SCALE WITH GAIN
rsp = prm.scale .* normrsp;                  % scale

pred = rsp;

%% COMPUTE ERROR

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end
end