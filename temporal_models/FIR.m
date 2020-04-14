function [err, pred] = FIR(param, data, stim, srate)
%
% function [err, pred] = FIR(param, data, stim, srate)
% INPUTS  -----------------------------------------------------------------
% params : 300 fields, reflecting the to-be-fitted linear filter
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

%% CONVOLVE THE IRF WITH THE STIMULUS
   
irf = param;
if size(param,2) ~=1, irf = irf'; end

linrsp  = conv2(stim, irf, 'full');         
linrsp  = linrsp(1:numtimepts,:);           

pred    = linrsp;

%% COMPUTE ERROR

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end
end