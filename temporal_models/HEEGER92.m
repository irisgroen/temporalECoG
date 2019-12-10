function [err, pred] = HEEGER92(param, data, stim, srate)

% Implementation of Heeger 1992 normalization model
% function [err, pred] = HEEGER92(param, data, stim, srate)
% INPUTS  -----------------------------------------------------------------
% params : 3 fields
%          1. tau1 -- time constant for linear impulse response
%          2. shift -- time between stimulus onset and when the signal reaches
%               the cortex (seconds)
%          3. sigma - semi-saturation constant
%          4. alpha - weight on feedback (0 = no feedback)
%          5. n     - output non-linearity exponent
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
fields = {'tau1', 'shift', 'sigma', 'alpha', 'n', 'rmax'};
prm      = toSetField([], fields, param);

%% Make impulse response function

t   = (1:numtimepts)' / srate;

irf = gammaPDF(t, prm.tau1, 2);

%% Initialize response

R = zeros(1, numtimepts); % Normalized response
G = zeros(size(R));       % Feedback signal
F = zeros(size(R));       % Multiplicative feedback
%K = 1;                   % Determines maximum responses
K = prm.rmax;               

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
    R(1) = max(L(1),0)^prm.n;
    G(1) = prm.alpha * R(1);
    F(1) = sqrt(K-G(1))/ prm.sigma;
    
    for ii = 2:length(t)
        R(ii) = max(L(ii) * F(ii-1),0)^prm.n;
        G(ii) = (1-prm.alpha) * G(ii-1) + prm.alpha * R(ii);
        G(ii) = min(G(ii), K);
        F(ii) = sqrt(K-G(ii-1)) / prm.sigma;
    end
    
    %pred(:,k) = R./max(R);
    pred(:,k) = R;
end

%% COMPUTE ERROR

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end

%%% FROM JINGS SIMULATION CODE (FIG 7 in plosCB paper)

% %% HEEGER 1992 MODEL
% 
% alpha = 0.01;
% delta = 3;
% R_max = 1; 
% sigma = 0.1;
% 
% R_t = zeros(size(t)); B_t = zeros(size(t));
% 
% for it = 2 : length(t)
%     B_t(it) = alpha * R_t(it) + (1- alpha) * B_t(max(it - delta, 1)); % this is where the exponential decay comes from
%     R_t(it + delta) = I(it)./sigma * (R_max - B_t(it));
% end
% 
% figure (1), subplot(7, 1, 3), cla, plot(t, normMax(R_t(1 : length(t))), 'r-', 'linewidth', 2.5), hold on, plot(t, normMax(I), 'k-'), axis tight, box off
% title('Heeger 1992'), set(gca, 'fontsize', 14)

end