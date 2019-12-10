function [params, pred] = tde_fitModel(objFunction, stim, data, srate, options, saveDir)

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
% <saveDir> directory to save parameters and fits, default
%
% NOTES
% TO DO add option to xvalidate on 'electrodes'? (train on all electrodes but 1, test on left out electrode) 
% This requires the function to know which electrodes belong to a xval set, so
% more inputs; or we should assume that the input always belongs to 1 set
 
%% PARSE OPTIONS

if ~exist('options', 'var') || isempty(options)
    options = struct();
end
if ~isfield(options,'startprm') || isempty(options.startprm)
	fprintf('[%s] Loading default starting values and bounds for model %s \n', mfilename, func2str(objFunction));
    %options.startprm = loadParamsFromJson(objFunction);
    options.startprm = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', func2str(objFunction))));
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
if ~exist('saveDir', 'var') || isempty(saveDir)
    saveDir = [];
end

% Model start point and bounds
x0 = options.startprm.x0;
lb = options.startprm.lb;
ub = options.startprm.ub;

% Check if plausible bounds are defined
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

% Set optimization options
searchopts = optimset('Display',options.display);
%searchopts.MaxFunEvals = options.maxiter;

%% FIT THE temporal model

% Initialize
nTimepts    = size(data,1);
nStim       = size(data,2);
nDatasets   = size(data,3);
nParams     = size(x0,2);
pred        = nan(nTimepts, nStim, nDatasets);
params      = nan(nParams,nDatasets);

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
    prm = nan(nParams,nFolds);
    prd = nan(nTimepts,nStim);
    
    for jj = 1:nFolds
        
        fit_inx = foldIndices(jj,:);
        pred_inx = setdiff(1:nFolds, fit_inx);
        if isempty(pred_inx), pred_inx = fit_inx; end
        
        stim2fit = stim(:,fit_inx);
        data2fit = dset(:,fit_inx);
        stim2predict = stim(:,pred_inx);
        
        fprintf('[%s] Fold %d: Fitting model on stimulus %s \n', mfilename, jj, num2str(fit_inx)) 
        fprintf('[%s] Fold %d: Predicting response for stimulus %s \n', mfilename, jj, num2str(pred_inx));
        
        % SEARCH FOR BEST-FITTING PARAMETERS
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
    params(:,ii) = mean(prm,2); % should this be mean or median? 

end

%% SAVE RESULTS
if ~isempty(saveDir)
    
    saveName = fullfile(saveDir, sprintf('%s_results', func2str(objFunction)));
    fprintf('[%s] Saving results to %s \n', mfilename, saveName);
      
    if exist(sprintf('%s.mat',saveName),'file')
        warning('[%s] Results file already exists! Writing new file with date-time stamp.',mfilename);
        saveName = sprintf('%s_%s', saveName, datestr(now,30));
        fprintf('[%s] Saving results to %s \n', mfilename, saveName);
    end
    save(saveName, 'pred', 'params');  
end

fprintf('[%s] Done!\n',mfilename);

end