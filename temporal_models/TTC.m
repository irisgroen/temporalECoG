function [err, pred] = TTC(param, data, stim, srate)
% Implementation of Two Temporal Channels model as implemented in Horiguchi
% et al, 2009, Neuroimage
% function [err, pred] = TTC(param, data, stim, srate)
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

%% USEFUL FUNCTIONS
normL2  = @(x) x./norm(x);
makeIRF = @(A, B, C, t)(t/A).^8 .* exp(-t/A) - 1 / B .* (t/C).^9 .* exp(-t/C);

%% SET UP THE MODEL PARAMETERS

% These values are from Table 1, row 1 (subject DT) in :
%   McKee and Taylor, JOSA 1984
% Presumably used by Horiguchi et al 2009
A = 3.29;
B = 14;
C = 3.85;

a = 2.75;
b = 11;
c = 3.18;

% sigma = 0.03; % smoothing kernel

fields = {'weight', 'shift', 'scale'};
prm      = toSetField([], fields, param);

%% COMPUTE THE IMPULSE RESPONSE FUNCTION

% make foveal channel impulse response
t   = 1000 * (1/srate : 1/srate : 0.100)'; % in milliseconds
irf_foveal = normL2(makeIRF(A, B, C, t));

% make peripheral channel impulse response
irf_peripheral = normL2(makeIRF(a, b, c, t));

% make transient and sustained IRFs as weighted sums of foveal and
%  peripheral IRFs
w = sum(irf_peripheral) / sum(irf_foveal );
irf_transient = normL2(irf_peripheral - w*irf_foveal);

w = 0.5;
irf_sustained = normL2(irf_foveal - w*irf_peripheral);

% debug: Compare to figure 7 (inset) in Horiguchi et al, 2009
% figure, plot(t, irf_sustained, 'b-', ...
%   t, irf_transient, 'k-', 'LineWidth', 3); xlim([0 100])

% % make smoothing kernel, transform from neuronal to ECoG broadband response
% s   = -1 : 1/srate : 1;
% ker = exp(-s.^2./(2 * sigma.^2));

%% COMPUTE THE NORMALIZATION RESPONSE

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

pred = rsp_combined;

%% COMPUTE ERROR

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end
end