function [err, pred] = dn_DNmodel(param, data, stim, srate)
%
% function [err, pred] = dn_DNmodel(param, data, stim, srate)
% INPUTS  -----------------------------------------------------------------
% params : 6 fields.
%          tau1 -- irf peak time, in unit of second
%          weight -- the weight in the biphasic irf function, set weight to
%          0 if want to use uniphasic irf function.
%          tau2 -- time window of adaptation, in unit of second
%          n -- exponent
%          sigma -- semi-saturation constant
%          shift -- time between stimulus onset and when the signal reaches
%          the cortex, in unit of second
%          scale -- response gain.
%
% data :   matrix, samples x trials
%
% stim :   matrix, samples x trial 
%
% srate :  samle rate in Hz
%
% OUTPUTS -----------------------------------------------------------------
% err:  sum of squared error
% pred: predicted time series


%% PRE-DEFINED /EXTRACTED VARIABLES

x       = []; % a struct of model parameteres

numtimepts  = size(stim,1);
numstim     = size(stim,2);

normSum = @(x) x./sum(x);

%% SET UP THE MODEL PARAMETERS

fields = {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};
x      = toSetField(x, fields, param);

%% COMPUTE THE IMPULSE RESPONSE FUNCTION

% HERE I ASSUME THAT THE NEGATIVE PART OF THE IMPULSE RESPONSE HAS A TIME
% CONSTANT 1.5 TIMES THAT OF THE POSITIVE PART OF THE IMPULSE RESPONSE
if x.tau1 > 0.5, warning('tau1>1, the estimation for other parameters may not be accurate'); end
    
t_irf   = 1/srate : 1/srate : 5;

irf_pos = gammaPDF(t_irf, x.tau1, 2);
irf_neg = gammaPDF(t_irf, x.tau1*1.5, 2);
irf     = irf_pos - x.weight.* irf_neg;

%% COMPUTE THE DELAYED REPSONSE FOR THE NORMALIZATION POOL

irf_norm = normSum(exp(-t_irf/x.tau2));

%% COMPUTE THE NORMALIZATION RESPONSE

linrsp  = NaN(size(stim));
numrsp  = NaN(size(stim));
poolrsp = NaN(size(stim));
demrsp  = NaN(size(stim));
normrsp = NaN(size(stim));

for istim = 1 : numstim
    % ADD SHIFT TO THE STIMULUS -------------------------------------------
    sft       = round(x.shift * srate);
    stimtmp   = padarray(stim(:, istim), [sft, 0], 0, 'pre');
    stim(:, istim) = stimtmp(1 : size(stim, 1));
    
    % COMPUTE THE NORMALIZATION NUMERATOR ---------------------------------
    linrsp(:, istim)  = convCut(stim(:, istim), irf, numtimepts);
    numrsp(:, istim)  = abs(linrsp(:, istim)).^x.n;
    
    % COMPUTE THE NORMALIZATION DENOMINATOR -------------------------------
    poolrsp(:, istim) = convCut(linrsp(:, istim), irf_norm, numtimepts);
    demrsp(:, istim)  = x.sigma.^x.n + abs(poolrsp(:, istim)).^x.n;
    
    % COMPUTE THE NORMALIZATION RESPONSE
    normrsp(:, istim) = x.scale.*(numrsp(:, istim)./demrsp(:, istim));
end

pred = normrsp;

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end
end