function [params, pred] = tde_fitModel(objFunction, stim, data, srate, options)

% function [params, pred] = tde_fitModel(objFunction, data, stim, srate, options) 
%
% <objFunction> model form
% <stimuli> stimulus time courses (time x condition)
% <data> provides the data as a cell vector of voxels x time.  can also be
%   X x Y x Z x time.  the number of time points should match the number of 
%   time points in <stimulus>.
% <srate> sample rate of the data
% <options> (optional) is a struct with the following fields:
%   <startprm>  starting values for and bounds for parameters (if not
%           provided, defaults are read in from a json file)
%   <algorithm>  search algorithm, either 'fminsearch' or 'bads' (default)
%   <xvalmode> method for cross validation (string), options are
%     'stimuli' : train on all stimulus conditions but 1, test on the left
%                   out stimulus
%     'none' : (default)
%   <display> is 'iter' | 'final' | 'off'.  default: 'iter'.
%
% NOTES
% TO DO add option to xvalidate on 'electrodes': train on all electrodes but 1, test on left out electrode? 
% requires function to know which electrodes belong to a group, so more inputs
 
%% FIT THE temporal model

% <options>
if ~exist('options', 'var') || isempty(options)
    options = struct();
end
if ~isfield(options,'startprm') || isempty(options.startprm)
	fprintf('[%s] Loading default starting values and bounds from json \n', mfilename);
    options.startprm = loadParamsFromJson(objFunction);
end
if ~isfield(options,'algorithm') || isempty(options.algorithm)
    options.algorithm = 'bads';
end
if ~isfield(options,'xvalmode') || isempty(options.xvalmode)
    options.xvalmode  = 'none';
end
if ~isfield(options,'display') || isempty(options.display)
    options.display   = 'iter';
end

% model start point and bounds
x0 = options.startprm.x0;
lb = options.startprm.lb;
ub = options.startprm.ub;

% check if plausible bounds are defined
switch options.algorithm
    case 'bads'
        if isfield(options.startprm, 'pub') 
            plb = options.startprm.plb;
            pub = options.startprm.pub;
        else
            fprintf('[%s] No plausible bounds specified for bads, switching to fminsearch \n', mfilename);
            options.algorithm = 'fminsearch';
        end
end

% initialize
nStim       = size(data,2);
nDatasets   = size(data,3);
pred        = nan(size(data));
params      = nan(size(x0,2),nDatasets);

searchopts = optimset('Display',options.display);
%searchopts.MaxFunEvals = options.maxiter;

for ii = 1:nDatasets % loop over channels or channel averages

    fprintf('[%s] Fitting model for dataset %d \n',mfilename, ii);
    
    dset = data(:,:,ii);

    %% SET UP CROSS VALIDATION
    switch options.xvalmode
       case 'stimuli'
           nFolds = nStim; % number of stimulus conditions
           foldIndices = nan(nStim, nStim-1);
           for jj = 1:nFolds 
               foldIndices(jj,:) = setdiff(1:nFolds,jj);
           end
       case 'none'
            nFolds = 1; foldIndices = 1:nStim;
        otherwise
            fprintf('[%s] Do not recognize this xvalmode, exiting \n',mfilename); return
    end
       
    %% FIT MODEL
    fprintf('[%s] Defined %d folds for fitting \n',mfilename, nFolds);
    prm = nan(size(x0,2),nFolds);
    prd = nan(size(dset));
    
    for jj = 1:nFolds
        
        fit_inx = foldIndices(jj,:);
        pred_inx = setdiff(1:nFolds, fit_inx);
        fprintf('[%s] Fitting model for fold %d on stimuli %s \n',mfilename, jj, num2str(fit_inx));

        if isempty(pred_inx), pred_inx = fit_inx; end
        
        stim2fit = stim(:,fit_inx);
        data2fit = dset(:,fit_inx);
        stim2predict = stim(:,pred_inx);
        
        switch options.algorithm
            case 'bads'
                prm(:,jj) = bads(@(x) objFunction(x, data2fit, stim2fit, srate),  x0, lb, ub, plb, pub, [], searchopts);
            case 'fminsearch'
                prm(:,jj) = fminsearchbnd(@(x) objFunction(x, data2fit, stim2fit, srate), x0, lb, ub, searchopts);
        end
        
        % GENERATE MODEL PREDICTION
        [~, prd(:,pred_inx)] = objFunction(prm(:,jj), [], stim2predict, srate);      
    end
    
    %% COLLECT FITTED PARAMETERS AND PREDICTIONS
    
    pred(:,:,ii) = prd;
    params(:,ii) = mean(prm,2);

end
end

%%% SUBROUTINES %%%

function options = loadParamsFromJson(objFunction)
    modelName = func2str(objFunction);
    default_opts = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', modelName)));
    options.x0 = default_opts.x0;
    options.lb = default_opts.lb;
    options.ub = default_opts.ub;
    if isfield(default_opts, 'pub')
        options.plb = default_opts.plb;
        options.pub = default_opts.pub;
    end
end