function [err, pred, numrsp, demrsp] = DN_modifiedIRF(param, data, stim, srate)
%
% function [err, pred, numrsp, demrsp] = DN(param, data, stim, srate)
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
% numrsp: numerator time course
% demrsp: denominator time course
%
% NOTE: the only difference with DN.m is that here the IRF for the
% numerator is also normalized to summed to 1, after the positive and
% negative IRF are combined (which are themselves each normalized to 1),
% and that the denominator and numerator are provided as output variables.

%% PRE-DEFINED /EXTRACTED VARIABLES
numtimepts  = size(stim,1);

%% USEFUL FUNCTIONS
normSum = @(x) x./sum(x);

%% SET UP THE MODEL PARAMETERS
fields = {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};
prm      = toSetField([], fields, param);

%% COMPUTE THE IMPULSE RESPONSE FUNCTION

% HERE I ASSUME THAT THE NEGATIVE PART OF THE IMPULSE RESPONSE HAS A TIME
% CONSTANT 1.5 TIMES THAT OF THE POSITIVE PART OF THE IMPULSE RESPONSE
    
t       = (1:numtimepts)' / srate;

irf_pos = gammaPDF(t, prm.tau1, 2);
irf_neg = gammaPDF(t, prm.tau1*1.5, 2);
irf     = irf_pos - prm.weight.* irf_neg;

% %%%% MODIFICATION: also sum numerator IRF to one
irf     = normSum(irf);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% COMPUTE THE DELAYED RESPONSE FOR THE NORMALIZATION POOL

irf_norm = normSum(exp(-t/prm.tau2));

%% COMPUTE THE NORMALIZATION RESPONSE

% ADD SHIFT TO THE STIMULUS -------------------------------------------
t2 = t-prm.shift;
stim = interp1(t, stim, t2, [], 0);

% COMPUTE THE NORMALIZATION NUMERATOR ---------------------------------
linrsp  = conv2(stim, irf, 'full');         % convolve
linrsp  = linrsp(1:numtimepts,:);           % cut
numrsp  = abs(linrsp).^prm.n;               % exponentiate

% COMPUTE THE NORMALIZATION DENOMINATOR -------------------------------
poolrsp = conv2(linrsp, irf_norm, 'full');  % convolve
poolrsp = poolrsp(1:numtimepts,:);          % cut
demrsp  = prm.sigma.^prm.n + abs(poolrsp).^prm.n; % exponentiate

% COMPUTE THE NORMALIZATION RESPONSE
normrsp = prm.scale.*(numrsp./demrsp);        % divide and scale

pred = normrsp;


%% COMPUTE ERROR

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end
end