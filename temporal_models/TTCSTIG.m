function [err, pred] = TTCSTIG(param, data, stim, srate)
%
% function [err, pred] = TTCSTIG(param, data, stim, srate)
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
x       = []; % a struct of model parameteres

numtimepts  = size(stim,1);
numstim     = size(stim,2);

normL2  = @(x) x./norm(x);

h = @(tau, n, t) (tau*factorial(n-1))^-1 * (t/tau).^(n-1) .* exp(-t/tau);

%% SET UP THE MODEL PARAMETERS


fields = {'weight', 'shift', 'scale'};
x      = toSetField(x, fields, param);


%% COMPUTE THE IMPULSE RESPONSE FUNCTION
t   = 1000 * (1/srate : 1/srate : 0.150)'; % in milliseconds

% make sustained channel impulse response
tau = 4.94;
n = 9;
irf_sustained = h(tau, n, t);

kappa = 1.33;
tau2  = kappa*tau;
n2    = 10;
xi    = 1.44;
irf_transient = xi*(irf_sustained - h(tau2, n2, t));



% debug: Compare to figure 7 (inset) in Horiguchi et al, 2009
% figure, plot(t_irf, irf_sustained, 'b-', ...
%   t_irf, irf_transient, 'k-', 'LineWidth', 3); xlim([0 .1])


% % make smoothing kernel, transform from neuronal to ECoG broadband response
% s   = -1 : 1/srate : 1;
% ker = exp(-s.^2./(2 * sigma.^2));

%% COMPUTE THE NORMALIZATION RESPONSE

rsp_transient  = NaN(size(stim));
rsp_sustained  = NaN(size(stim));
rsp_combined = NaN(size(stim));

for istim = 1 : numstim
    % ADD SHIFT TO THE STIMULUS -------------------------------------------
    sft       = round(x.shift * srate);
    stimtmp   = padarray(stim(:, istim), [sft, 0], 0, 'pre');
    stim(:, istim) = stimtmp(1 : size(stim, 1));
    
    % COMPUTE THE TRANSIENT RESPONSE ---------------------------------
    rsp_transient(:, istim)  = convCut(stim(:, istim), irf_transient, numtimepts);
    rsp_transient(:, istim)  = abs(rsp_transient(:, istim)).^2;
    
    % COMPUTE THE SUSTAINED RESPONSE -------------------------------
    rsp_sustained(:, istim) = convCut(stim(:, istim), irf_sustained, numtimepts);
    
    % COMPUTE THE COMBINED RESPONSE
    rsp_combined(:, istim) = x.scale.*(x.weight*rsp_transient(:, istim) ...
        + (1-x.weight)*rsp_sustained(:, istim));
end

pred = rsp_combined;

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end
end