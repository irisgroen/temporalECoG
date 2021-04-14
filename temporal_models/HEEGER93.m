function [err, pred] = HEEGER93(param, data, stim, srate)

% Implementation of Heeger 1993 normalization model
% function [err, pred] = HEEGER93(param, data, stim, srate)
% INPUTS  -----------------------------------------------------------------
% params : 5 fields
%          1. tau1  -- time constant for linear impulse response
%          2. shift -- time between stimulus onset and when the signal reaches
%               the cortex (seconds)
%          3. sigma -- semi-saturation constant
%          4. alpha -- weight on feedback (0 = no feedback)
%          5. rmax  -- maximum response
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
%
% NOTE: the only difference with the HEEGER92 model is the exponent in the
% denominator (2 instead of 1).


%% PRE-DEFINED /EXTRACTED VARIABLES
numtimepts  = size(stim,1);
numstim     = size(stim,2);
n           = 2;  % exponent               

%% SET UP THE MODEL PARAMETERS
fields = {'tau1', 'shift', 'sigma', 'alpha',  'rmax'};
prm      = toSetField([], fields, param);

%% Make impulse response function

t   = (1:numtimepts)' / srate;
irf = gammaPDF(t, prm.tau1, 2);

%% Initialize response

R = zeros(numtimepts, numstim); % Normalized response
G = zeros(size(R));             % Feedback signal
F = zeros(size(R));             % Multiplicative feedback
K = prm.rmax;                   % Determines maximum response

%% Compute normalization response

% ADD SHIFT TO THE STIMULUS -------------------------------------------
t2 = t-prm.shift;
stim = interp1(t, stim, t2, [], 0);

% COMPUTE THE NORMALIZATION NUMERATOR ---------------------------------
linrsp  = conv2(stim, irf, 'full');         % convolve
linrsp  = linrsp(1:numtimepts,:);           % cut

% COMPUTE RESPONSE AT EACH TIMEPOINT 
R(1,:) = max(linrsp(1,:),0).^n;
G(1,:) = prm.alpha * R(1,:);
F(1,:) = sqrt(K-G(1,:)) / prm.sigma;

for ii = 2:length(t)
    R(ii,:) = max(linrsp(ii,:) .* F(ii-1,:),0).^n;
    G(ii,:) = (1-prm.alpha) * G(ii-1,:) + prm.alpha * R(ii,:);
    G(ii,:) = min(G(ii,:), K);
    F(ii,:) = sqrt(K-G(ii-1,:)) / prm.sigma;
end

pred = R;

%% COMPUTE ERROR

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end

end