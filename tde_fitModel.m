function [params, pred] = tde_fitModel(objFunction, data, stim, srate, opts)
% Description
% To call fitting function, we need:
%   1. the objective function (model form and type of error)
%   2. data
%   3. stimuli
%   4. starting values and bounds for parameters (there should be defaults
%               for each type of model)
%   5. sample rate of the data
 
%% FIT THE temporal model

% <opts>
if ~exist('opts', 'var') || isempty(opts)
	sprintf('Loading default opts from json');
    modelName = func2str(objFunction);
    opts = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', modelName)));
end 

% model start point and bounds
x0 = opts.x0;
lb = opts.lb;
ub = opts.ub;

% check if plausible bounds are defined
useBads = 0;
if isfield(opts, 'pub')
    plb = opts.plb;
    pub = opts.pub;
    useBads = 1;
end

% initialize
nDatasets   = size(data,3);
pred        = nan(size(data));
params      = nan(size(opts.x0,2),nDatasets);

options = optimset('Display','iter');
% options.MaxFunEvals = 10000;

for ii = 1:nDatasets % loop over channels or channel averages

    fprintf('[%s] Fitting model for dataset %d \n',mfilename, ii);
    
    data2fit = data(:,:,ii);
         
    %% FIT MODEL

    % use bads if plausible upper and lower bounds are defined, otherwise use fminsearch
    
    if useBads
        prm = bads(@(x) objFunction(x, data2fit, stim, srate),  x0, lb, ub, plb, pub, [], options);
    else
        prm = fminsearchbnd(@(x) objFunction(x, data2fit, stim, srate), x0, lb, ub, options);
    end
    
    %% GENERATE MODEL PREDICTIONS

    [~, p] = objFunction(prm, [], stim, srate);
    % p = p./max(p(:));
    pred(:,:,ii) = p;
    
    %% COLLECT FITTED PARAMETERS
    
    params(:,ii) = prm;

end