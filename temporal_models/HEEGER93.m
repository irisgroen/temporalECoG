function [err, pred] = HEEGER93(param, data, stim, srate)

% Implementation of Heeger 1992 normalization model
% function [err, pred] = HEEGER92(param, data, stim, srate)
% INPUTS  -----------------------------------------------------------------
% params : 3 fields
%          1. tau1 -- time constant for linear impulse response
%          2. shift -- time between stimulus onset and when the signal reaches
%               the cortex (seconds)
%          3. sigma - semi-saturation constant
%          4. alpha - weight on feedback (0 = no feedback)
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
fields = {'tau1', 'shift', 'sigma', 'alpha',  'rmax'};
prm      = toSetField([], fields, param);

%% Make impulse response function

t   = (1:numtimepts)' / srate;
irf = gammaPDF(t, prm.tau1, 2);

%% Initialize response

R = zeros(1, numtimepts); % Normalized response
G = zeros(size(R));       % Feedback signal
F = zeros(size(R));       % Multiplicative feedback
K = prm.rmax;             % Determines maximum responses
               

%% Compute normalization response

% ADD SHIFT TO THE STIMULUS -------------------------------------------
sft       = round(prm.shift * srate);
stimtmp   = padarray(stim, [sft, 0], 0, 'pre');
stim      = stimtmp(1 : size(stim, 1), :);

% COMPUTE THE NORMALIZATION NUMERATOR ---------------------------------
linrsp  = conv2(stim, irf, 'full');         % convolve
linrsp  = linrsp(1:numtimepts,:);           % cut

pred = [];
for k = 1 : size(stim, 2)
    %L    = convCut(stim(k, :), irf, length(stim));
    L = linrsp(:,k);
    R(1) = max(L(1),0)^2;
    G(1) = prm.alpha * R(1);
    F(1) = sqrt(K-G(1))/ prm.sigma;
    
    for ii = 2:length(t)
        R(ii) = max(L(ii) * F(ii-1),0)^2;
        F(ii) = sqrt(K-G(ii-1)) / prm.sigma;
        G(ii) = (1-prm.alpha) * G(ii-1) + prm.alpha * R(ii);
        G(ii) = min(G(ii), K);
    end
    
    pred(:,k) = R;
end

%% COMPUTE ERROR

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end

end