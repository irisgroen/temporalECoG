function [params, pred] = tde_fitModel(objFunction, stim, data, srate, options, saveDir, saveName)

% [params, pred] = tde_fitModel(objFunction, data, stim, srate, options, saveDir, saveName) 
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
%     0 : no cross-validation (default)
%     1 : train on all stimulus conditions but 1, test on the left out 
%   <display> is 'iter' | 'final' | 'off'.  default: 'iter'.
% <saveDir> path to save parameters and fits; if empty, results are not
%   saved (default)
% <saveName> string to add to the save filename, if results are saved
%   (default empty)
%
% Example

%% PARSE OPTIONS

if ~exist('options', 'var') || isempty(options), options = struct(); end
if ~isfield(options,'algorithm') || isempty(options.algorithm), options.algorithm = 'bads'; end
if ~isfield(options,'xvalmode') || isempty(options.xvalmode), options.xvalmode = 0; end
if ~isfield(options,'display') || isempty(options.display), options.display = 'iter'; end
if ~exist('saveDir', 'var'), saveDir = []; end
if ~exist('saveName', 'var'), saveName = []; end
if iscell(objFunction), objFunction = objFunction{1}; end

% Get model start points and bounds
if ~isfield(options,'startprm') || isempty(options.startprm)
	fprintf('[%s] Loading default starting values and bounds for model %s \n', mfilename, func2str(objFunction));
    options.startprm = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', func2str(objFunction))));
end
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

%% SET UP CROSS VALIDATION
switch options.xvalmode
    case 0
        nFolds = 1; foldIndices{1} = 1:nStim;
        fprintf('[%s] Xval mode is 0: generating predictions by fitting the full dataset \n',mfilename);
    case 1
        nFolds = nStim+1; % number of stimulus conditions + full
        foldIndices = cell(nFolds,1); % first fold is full
        foldIndices{1} = 1:nStim;
        for jj = 1:nStim, foldIndices{1+jj} = setdiff(1:nStim,jj); end
        fprintf('[%s] Xval mode is 1: generating predictions for left-out-stimulus across %d folds \n',mfilename, nFolds-1);
    otherwise
        fprintf('[%s] Xval mode not recognized, exiting \n',mfilename); return
end

%% FIT MODEL

for ii = 1:nDatasets % loop over channels or channel averages

    dset = data(:,:,ii);
    
    fprintf('[%s] Fitting model for dataset %d \n',mfilename, ii);
       
    for jj = 1:nFolds
        
        fit_inx = foldIndices{jj};
        pred_inx = setdiff(1:nStim, fit_inx);
        if isempty(pred_inx), pred_inx = fit_inx; end
        
        stim2fit = stim(:,fit_inx);
        data2fit = dset(:,fit_inx);
        stim2predict = stim(:,pred_inx);
        
%         if jj > 1
%             fprintf('[%s] Fold %d: Fitting on stimulus %s \n', mfilename, jj-1, num2str(fit_inx)) 
%             fprintf('[%s] Fold %d: Predicting for stimulus %s \n', mfilename, jj-1, num2str(pred_inx));
%         end
        
        % Search for best-fitting parameters
        switch options.algorithm
            case 'bads'
                prm = bads(@(x) objFunction(x, data2fit, stim2fit, srate),  x0, lb, ub, plb, pub, [], searchopts);
            case 'fminsearch'
                prm = fminsearchbnd(@(x) objFunction(x, data2fit, stim2fit, srate), x0, lb, ub, searchopts);
        end
        
        % Save params from full model fit
        if jj == 1, params(:,ii) = prm; end
        
        % Generate model prediction
        [~, pred(:,pred_inx,ii)] = objFunction(prm, [], stim2predict, srate);      
    end
end

%% SAVE RESULTS
if ~isempty(saveDir)
    
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    if isempty(saveName)
        saveName = sprintf('%s_xvalmode%d', func2str(objFunction), options.xvalmode);
    else
        saveName = sprintf('%s_xvalmode%d_%s', func2str(objFunction), options.xvalmode, saveName);
    end
    saveName = fullfile(saveDir, saveName);
    fprintf('[%s] Saving results to %s \n', mfilename, saveName);
      
    if exist(sprintf('%s.mat',saveName),'file')
        warning('[%s] Results file already exists! Writing new file with date-time stamp.',mfilename);
        saveName = sprintf('%s_%s', saveName, datestr(now,30));
        fprintf('[%s] Saving results to %s \n', mfilename, saveName);
    end
    save(saveName, 'pred', 'params', 'stim', 'data', 'srate', 'options', 'objFunction');  
end

fprintf('[%s] Done!\n',mfilename);

end