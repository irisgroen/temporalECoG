function [results] = tde_evaluateModelFit(data, objFunction, params, pred, stim_info)

% Inputs: cell arrays with objFunction, params, data, pred
if ~iscell(objFunction), objFunction = {objFunction}; end
if ~iscell(params), params = {params}; end
if ~iscell(pred), pred = {pred}; end

nModels     = size(objFunction,2);
[~,nStim,nDatasets] = size(data);
stimcond    = stim_info.condition;

% initialize
R2stim      = nan(nStim,nDatasets);
R2concat    = nan(1,nDatasets);
R2cond      = nan(length(unique(stimcond)), nDatasets); 

derivedPrm  = nan(2,nDatasets);
derivedPred = [];

%% COMPUTE SUMMARY METRICS
results = struct;

% Loop over models
for kk = 1:nModels
    
    fprintf('[%s] Evaluating model %s for %d datasets\n', mfilename, func2str(objFunction{kk}), nDatasets)
    
    % Loop over channels or channel averages
    for ii = 1:nDatasets
        
        %% COMPUTE R-SQUARE
        
        % For each individual stimulus 
        rsq = nan(nStim,1);
        for jj = 1:nStim % loop over stimuli
            DATA = data(:,jj,ii);
            MODEL = pred{kk}(:,jj,ii);
            rsq(jj) = computeR2(DATA,MODEL);
        end
        R2stim(:,ii) = rsq;
             
        % For all stimuli concatenated
        DATA = flatten(data(:,:,ii));
        MODEL = flatten(pred{kk}(:,:,ii));
        R2concat(:,ii) = computeR2(DATA,MODEL);
        
        % For specific conditions
        cond = unique(stimcond); nCond = length(cond);
        for jj = 1:nCond
            condInx = find(stimcond == cond(jj));
            DATA = flatten(data(:,condInx,ii));
            MODEL = flatten(pred{kk}(:,condInx,ii));
            R2cond(jj,ii) = computeR2(DATA,MODEL);
        end
        
        fprintf('[%s] R2 for dataset %d = %0.2f \n', mfilename, ii, R2concat(:,ii))
    

        %% COMPUTE MODELBASED DERIVED PARAMETERS 
        fprintf('[%s] Computing derived parameters...\n', mfilename)

        % Generate a prediction to a sustained stimulus:
        [derived_prm, pred_derived] = tde_computeDerivedParams(objFunction{kk}, params{kk}(:,ii));

        derivedPrm(1,ii) = derived_prm.t2pk;
        derivedPrm(2,ii) = derived_prm.r_asymp;
        derivedPrm(3,ii) = derived_prm.r_double;
        derivedPrm(4,ii) = derived_prm.t_isi;

        derivedPred(:,ii)= pred_derived;
      
    end
    
    %% COLLECT RESULTS
    results(kk).model           = objFunction{kk};
	results(kk).params          = params{kk};
    results(kk).pred            = pred{kk};
    results(kk).R2.stim         = R2stim;
    results(kk).R2.concat_all   = R2concat;
    results(kk).R2.concat_cond  = R2cond; 
    results(kk).derived.names   = {'T2peak', 'Rasymptote', 'Rdouble', 'Tisi'};
    results(kk).derived.params  = derivedPrm;
    results(kk).derived.pred    = derivedPred;
end
fprintf('[%s] Done!\n', mfilename)
end

