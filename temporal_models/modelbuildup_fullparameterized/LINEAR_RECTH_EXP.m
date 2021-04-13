function [err, pred] = LINEAR_RECTH_EXP(param, data, stim, srate)
%
% function [err, pred] = LINEAR_RECTH_EXP(param, data, stim, srate)
% INPUTS  -----------------------------------------------------------------
% params : 7 fields.
%          1. tau1 -- time to peak for positive IRF (seconds)
%          2. n1 -- phase delay for positive IRF 
%          3. tau2 -- time to peak for negative IRF (seconds)
%          4. n2 -- phase delay for negative IRF 
%          5. weight -- ratio of negative to positive IRFs 
%          6. shift -- time between stimulus onset and when the signal reaches
%               the cortex (seconds)
%          7. scale -- response gain.
%          8. n -- exponent 
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
fields = {'tau1', 'n1', 'tau2', 'n2', 'weight', 'shift', 'scale', 'n'};
prm      = toSetField([], fields, param);

%% COMPUTE THE IMPULSE RESPONSE FUNCTION

t       = (1:numtimepts)' / srate;

% COMPUTE THE IRF
n1 = round(prm.n1); % n has to be an integer
n2 = round(prm.n2); % n has to be an integer
if n2 == 0, n2 = 1; end
if n1 == 0, n1 = 1; end
irf_pos = gammaPDF(t, prm.tau1, n1);
irf_neg = gammaPDF(t, prm.tau2, n2);
irf     = irf_pos - prm.weight.* irf_neg;

%% COMPUTE THE RESPONSE

% ADD SHIFT TO THE STIMULUS -------------------------------------------
t2 = t-prm.shift;
stim = interp1(t, stim, t2, [], 0);

% COMPUTE THE CONVOLUTION
linrsp  = conv2(stim, irf, 'full');         % convolve
linrsp  = linrsp(1:numtimepts,:);           % cut
linrsp  = max(linrsp,0); % half wave rectification
linrsp  = linrsp.^prm.n; % exponentiate

% SCALE WITH GAIN
rsp = prm.scale .* linrsp;                    % scale

pred = rsp;

%% COMPUTE ERROR

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end
end