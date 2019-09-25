function [err, pred] = dn_DNmodel(param, data, stim, t)
%
% function [err, pred] = dn_DNmodel(param, data, stim, t)
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
% OUTPUTS -----------------------------------------------------------------
% err: sum of squared error.
% pred: predicted time series


%% PRE-DEFINED /EXTRACTED VARIABLES

x       = []; % a struct of model parameteres

t_lth   = length(t);
dt      = t(2) - t(1);

normSum = @(x) x./sum(x);

%% SET UP THE MODEL PARAMETERS

fields = {'tau1', 'weight', 'tau2', 'n', 'sigma', 'shift', 'scale'};
x      = toSetField(x, fields, param);

%% COMPUTE THE IMPULSE RESPONSE FUNCTION

% HERE I ASSUME THAT THE NEGATIVE PART OF THE IMPULSE RESPONSE HAS A TIME
% CONSTANT 1.5 TIMES THAT OF THE POSITIVE PART OF THE IMPULSE RESPONSE
if x.tau1 > 0.5, warning('tau1>1, the estimation for other parameters may not be accurate'); end
    
t_irf   = dt : dt : 5;

irf_pos = gammaPDF(t_irf, x.tau1, 2);
irf_neg = gammaPDF(t_irf, x.tau1*1.5, 2);
irf     = irf_pos - x.weight.* irf_neg;

%% COMPUTE THE DELAYED REPSONSE FOR THE NORMALIZATION POOL

irf_norm = normSum(exp(-t_irf/x.tau2));

%% COMPUTE THE NORMALIZATION RESPONSE

for istim = 1 : size(stim, 1)
    % ADD SHIFT TO THE STIMULUS -------------------------------------------
    sft       = round(x.shift * (1/dt));
    stimtmp   = padarray(stim(istim, :), [0, sft], 0, 'pre');
    stim(istim, :) = stimtmp(1 : size(stim, 2));
    
    % COMPUTE THE NORMALIZATION NUMERATOR ---------------------------------
    linrsp(istim, :)  = convCut(stim(istim, :), irf, t_lth);
    numrsp(istim, :)  = linrsp(istim, :).^x.n;
    
    % COMPUTE THE NORMALIZATION DENOMINATOR -------------------------------
    poolrsp(istim, :) = convCut(linrsp(istim, :), irf_norm, t_lth);
    demrsp(istim, :)  = x.sigma.^x.n + poolrsp(istim, :).^x.n;
    
    % COMPUTE THE NORMALIZATION RESPONSE
    normrsp(istim, :) = x.scale.*(numrsp(istim, :)./demrsp(istim, :));
end

pred = normrsp;

if isempty(data)
    err = []; 
else
    err = sum((pred(:) - data(:)).^2);
end
end