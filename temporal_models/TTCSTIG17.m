function [err, pred] = TTCSTIG17(param, data, stim, srate)
%
% function [err, pred] = TTCSTIG17(param, data, stim, srate)
% INPUTS  -----------------------------------------------------------------
% params : 3 fields
%          1. weight - relative weight on transient channel ([0 1]) 
%          2. shift -- time between stimulus onset and when the signal reaches
%               the cortex (seconds)
%          3. scale -- response gain.
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
t_stim      = (1:numtimepts)' / srate;

%% FUNCTIONS
h = @(tau, n, t) (tau*factorial(n-1))^-1 * (t/tau).^(n-1) .* exp(-t/tau);

%% SET UP THE MODEL PARAMETERS
fields = {'weight', 'shift', 'scale'};
prm      = toSetField([], fields, param);

%% COMPUTE THE IMPULSE RESPONSE FUNCTION
t   = 1000 * (1/srate : 1/srate : 0.150)'; % in milliseconds

% make sustained channel impulse response
tau = 4.94;
n1 = 9;
irf_sustained = h(tau, n1, t);

% make transient channel impulse response
kappa = 1.33;
tau2  = kappa*tau;
n2    = 10;
xi    = 1.44;
irf_transient = xi*(irf_sustained - h(tau2, n2, t));

% debug: 
% figure, plot(t_irf, irf_sustained, 'b-', ...
%   t_irf, irf_transient, 'k-', 'LineWidth', 3); xlim([0 .1])

%% COMPUTE THE PREDICTED TWO CHANNEL RESPONSE

% ADD SHIFT TO THE STIMULUS -------------------------------------------
t2 = t_stim-prm.shift;
stim = interp1(t_stim, stim, t2, [], 0);
  
% COMPUTE THE TRANSIENT RESPONSE --------------------------------------
rsp_transient = conv2(stim, irf_transient, 'full'); % convolve
rsp_transient = rsp_transient(1:numtimepts,:);      % cut
rsp_transient = abs(rsp_transient).^2;              % square

% COMPUTE THE SUSTAINED RESPONSE -------------------------------
rsp_sustained = conv2(stim, irf_sustained, 'full'); % convolve
rsp_sustained = rsp_sustained(1:numtimepts,:);      % cut
      
% COMPUTE THE COMBINED RESPONSE
rsp_combined = prm.scale.*(prm.weight*rsp_transient ...
        + (1-prm.weight)*rsp_sustained);

%% COMPUTE ERROR

pred = rsp_combined;

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end
end