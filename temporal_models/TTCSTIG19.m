function [err, pred] = TTCSTIG19(param, data, stim, srate)
%
% function [err, pred] = TTCSTIG19(param, data, stim, srate)
% INPUTS  -----------------------------------------------------------------
% params : 3 fields
%          1. weight -- relative weight on transient channel ([0 1]) 
%          2. shift -- time between stimulus onset and when the signal reaches
%               the cortex (seconds)
%          3. scale -- response gain.
%          4. tau - time constant for neural IRFs
%          5. alpha -- adaptation controlling exponential decay of sustained channel
%          6. k_on -- sigmoid on nonlinearity of transient channel
%          7. k_off -- sigmoid off nonlinearity of transient channel
%          8. lambda - scale of transient channel 
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


%%% FROM TCN toolbox
% % sustained response: (stimulus * sustained IRF) x exponential[tau_ae]
% conv_snS = @(s, tau_s, tau_ae) cellfun(@(X, Y, ON, OFF) code_exp_decay(X, ON, OFF, Y, fs), ...
%     cellfun(@(XX, YY) convolve_vecs(XX, YY, 1, 1), s, repmat({nrfS_fun(tau_s)}, nruns, 1), 'uni', false), ...
%     repmat({adapt_fun(tau_ae)}, nruns, 1), model.onsets, model.offsets, 'uni', false);
% % transient response: tch_sigmoid(stimulus * transient IRF)
% conv_snT = @(s, tau_s) cellfun(@(X, Y) convolve_vecs(X, Y, 1, 1), ...
%     s, repmat({nrfT_fun(tau_s)}, nruns, 1), 'uni', false);
% conv_snTs = @(s, tau_s, Lp, Kp, Kn) cellfun(@(X, lp, kp, kn) tch_sigmoid(X, lp, kp, lp, kn), ...
%     conv_snT(s, tau_s), repmat({Lp}, nruns, 1), repmat({Kp}, nruns, 1), repmat({Kn}, nruns, 1), 'uni', false);
%%%

%% PRE-DEFINED /EXTRACTED VARIABLES
numtimepts  = size(stim,1);

%% USEFUL FUNCTIONS
normL2  = @(x) x./norm(x);
h = @(tau, n, t) (tau*factorial(n-1))^-1 * (t/tau).^(n-1) .* exp(-t/tau);

%% SET UP THE MODEL PARAMETERS
fields = {'weight', 'shift', 'scale', 'tau', 'k_on', 'k_off', 'lambda', 'tau_ae'};
prm      = toSetField([], fields, param);

%% COMPUTE THE IMPULSE RESPONSE FUNCTION
t   = 1000 * (1/srate : 1/srate : 0.150)'; % in milliseconds

% make sustained channel impulse response
tau = prm.tau; % tau = 4.94; % fitted here instead of fixed in STIG17
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
sft       = round(prm.shift * srate);
stimtmp   = padarray(stim, [sft, 0], 0, 'pre');
stim = stimtmp(1 : size(stim, 1), :);
  
% COMPUTE THE TRANSIENT RESPONSE --------------------------------------
rsp_transient = conv2(stim, irf_transient, 'full'); % convolve
rsp_transient = rsp_transient(1:numtimepts,:);      % cut
%rsp_transient = abs(rsp_transient).^2;              % square
% instead, multiply with sigmoid:
rsp_transient = tch_sigmoid(rsp_transient, prm.lambda, prm.k_on, prm.lambda, prm.k_off);

% COMPUTE THE SUSTAINED RESPONSE -------------------------------
rsp_sustained = conv2(stim, irf_sustained, 'full'); % convolve
rsp_sustained = rsp_sustained(1:numtimepts,:);      % cut

% multiply with adaptation function
adapt_exp = exp(-(1:60000) / prm.tau_ae);
starts = 0; stops = numtimepts *(1/srate);% should be updated with actual start/stop times of stim
rsp_sustained = code_exp_decay(rsp_sustained, starts, stops, adapt_exp, srate);

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

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = tch_sigmoid(X, lambda_p, kappa_p, lambda_n, kappa_n)

if nargin < 4 || isempty(lambda_n); lambda_n = lambda_p; end
if nargin < 5 || isempty(kappa_n); kappa_n = kappa_p; end
X = rectify(X, 'abs', .001);
X_p = X; X_p(X < 0) = 0;
X_n = X; X_n(X_n > 0) = 0;
weibull_p = 1 - exp(-(X_p ./ lambda_p) .^ kappa_p);
weibull_n = 1 - exp(-(-X_n ./ lambda_n) .^ kappa_n);
Y = weibull_p + weibull_n;

end

function resp_decay = code_exp_decay(resp_in, starts, stops, decay_exp, fs)

    % Helper function for coding stimulus-specific exponential response decay. 
%
% INPUTS
%   1) resp_in: input activity matrix (frames x predictors)
%   2) starts: beginnings of decay activity windows (seconds)
%   3) stops: ends of decay activity windows (seconds)
%   4) decay_exp: exponential function modeling decay of activity
%   5) fs: temporal sampling rate of resp_in and resp_out (Hz)
%
% OUTPUT
%   resp_decay: output activity matrix with decay (frames x predictors)
%
% AS 10/2017

if ~isempty(resp_in)
    decay_fun = zeros(size(resp_in, 1), 1);
    for ss = 1:length(starts)
        start_idx = round(starts(ss) * fs); stop_idx = round(stops(ss) * fs);
        decay_idxs = start_idx + 1:stop_idx; dl = length(decay_idxs);
        decay_idxs = decay_idxs(1:min([dl length(decay_exp)]));
        decay_fun(decay_idxs) = decay_exp(1:min([dl length(decay_exp)]));
    end
    resp_decay = resp_in .* repmat(decay_fun, 1, size(resp_in, 2));
else
    resp_decay = [];
end


end
end