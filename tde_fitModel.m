function [params, pred, pnames] = tde_fitModel(objFunction, stim, data, srate, options)

% [params, pred, pnames] = tde_fitModel(objFunction, data, stim, srate, options) 
%
% Fits a temporal model to neural response time courses
%
% <objFunction> model form
% <stimuli> stimulus time courses (time x condition)
% <data> the data as a vector of time x conditions x electrodes. The number
%   of time points should match the number of time points in <stim>.
% <srate> sample rate of the data
% <options> (optional) is a struct with the following fields:
%   <startprm>  starting values for and bounds for parameters (if not
%           provided, defaults are read in from a json file)
%   <algorithm>  search algorithm, options are 
%     'fmincon' 
%     'lsqnonlin' 
%     'bads' (default)
%     'surrogateopts' 
%   note that fmincon and lsqnonlin can't use integerconstraints (necessary
%   when fitting the n going into the factoral of the IRF, see gammaPDF.m)
%   <xvalmode> method for cross validation (string), options are
%     0 : no cross-validation (default)
%     1 : train on all stimulus conditions but 1, test on the left out 
%   <display> is 'iter' | 'final' | 'off'.  default: 'iter'.
%
% 2022 Iris Groen

%% PARSE OPTIONS

if ~exist('options', 'var') || isempty(options), options = struct(); end
% Model fitting options
if ~isfield(options,'algorithm') || isempty(options.algorithm), options.algorithm = 'bads'; end
if ~isfield(options,'xvalmode') || isempty(options.xvalmode), options.xvalmode = 0; end
if ~isfield(options,'display') || isempty(options.display), options.display = 'iter'; end

% Some formatting
if iscell(objFunction), objFunction = objFunction{1}; end

% Get model start points and bounds
if ~isfield(options,'startprm') || isempty(options.startprm)
	fprintf('[%s] Loading default starting values and bounds for model %s \n', mfilename, func2str(objFunction));
    options.startprm = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', func2str(objFunction))));
end
x0 = options.startprm.x0;
lb = options.startprm.lb;
ub = options.startprm.ub;
pnames = strsplit(options.startprm.params,',');

% Check if plausible bounds are defined
switch options.algorithm
    case 'bads'
        if isfield(options.startprm, 'pub') 
            plb = options.startprm.plb;
            pub = options.startprm.pub;
        else
            fprintf('[%s] No plausible bounds specified for bads, switching to lsqnonlin \n', mfilename);
            options.algorithm = 'lsqnonlin';
        end
end

fprintf('[%s] Algorithm = %s \n', mfilename, options.algorithm);

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
    
    fprintf('[%s] Fitting model for dataset %d out of %d \n',mfilename, ii, nDatasets);
       
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
                % Set optimization options
                searchopts = optimset('Display',options.display);
                searchopts.MaxIterations = 10000;
                searchopts.MaxFunctionEvaluations = 10000;
                prm = bads(@(x) objFunction(x, data2fit, stim2fit, srate),  x0, lb, ub, plb, pub, [], searchopts);
            case 'lsqnonlin'
                searchopts.Algorithm = 'levenberg-marquardt';
                searchopts = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
                searchopts.MaxIterations = 10000;
                searchopts.MaxFunctionEvaluations = 10000;
                searchopts.Display = options.display;
                prm = lsqnonlin(@(x) objFunction(x, data2fit, stim2fit, srate), x0, lb, ub, searchopts);
            case 'fmincon'
                searchopts = optimoptions(@fmincon, 'Algorithm', 'sqp');
                searchopts.MaxIterations = 10000;
                searchopts.MaxFunctionEvaluations = 10000;
                searchopts.Display = options.display;
                problem.objective = @(x) objFunction(x, data2fit, stim2fit, srate);
                problem.x0 = x0;
                problem.solver = 'fmincon';
                problem.lb = lb;
                problem.ub = ub;
                problem.options = searchopts;
                prm = fmincon(problem);
            case 'surrogateopt' 
                searchopts = optimoptions('surrogateopt','PlotFcn',[], "ConstraintTolerance",1e-6);
                searchopts.MaxFunctionEvaluations = 10000;
                searchopts.Display = options.display;
                intcon = contains(pnames, 'n_irf');
                if any(intcon), intcon_idx = find(intcon); else, intcon_idx = []; end
                prm = surrogateopt(@(x) objFunction(x, data2fit, stim2fit, srate),lb,ub,intcon_idx,searchopts);       
            otherwise
                error('[%s] Fitting algorithm not recognized \n', mfilename);
        end
        
        % Save params from full model fit
        if jj == 1, params(:,ii) = prm; end
        
        % Generate model prediction
        [~, pred(:,pred_inx,ii)] = objFunction(prm, [], stim2predict, srate);      
    end
end

fprintf('[%s] Done with model fitting!\n',mfilename);

end