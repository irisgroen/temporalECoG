function [results] = tde_evaluateModelFit(d, includeDerivedParams)
 
% Input: struct with data and fitted model parameters and predictions

if ~exist('includeDerivedParams', 'var') || isempty(includeDerivedParams)
    includeDerivedParams = true;
end

nModels = size(d,2);

%% COMPUTE SUMMARY METRICS
results = struct;

% Loop over models
for kk = 1:nModels
    
    data        = d(kk).data;
    pred        = d(kk).pred;
    params      = d(kk).params;
    objFunction = d(kk).objFunction;
    stim_info   = d(kk).stim_info;
    stimcond    = stim_info.condition;
    
    % Get sizes of data and stiminfo
    [~,nStim,nDatasets] = size(data);

    fprintf('[%s] Evaluating model %s for %d datasets\n', mfilename, func2str(objFunction), nDatasets)

    % Initialize
    R2stim       = nan(nStim,nDatasets);
    R2concat     = nan(1,nDatasets);
    R2cond       = nan(length(unique(stimcond)), nDatasets); 
    derivedPrm   = nan(6,nDatasets);
    derivedPredS = nan(2000,nDatasets);
    derivedPredT = nan(500,nDatasets);
    pred_names   = [];
    
    % Loop over channels or channel averages
    for ii = 1:nDatasets
        
        %% COMPUTE R-SQUARE
        
        % For each individual stimulus 
        rsq = nan(nStim,1);
        for jj = 1:nStim % loop over stimuli
            DATA = data(:,jj,ii);
            MODEL = pred(:,jj,ii);
            rsq(jj) = computeR2(DATA,MODEL);
        end
        R2stim(:,ii) = rsq;
             
        % For all stimuli concatenated
        DATA = flatten(data(:,:,ii));
        MODEL = flatten(pred(:,:,ii));
        R2concat(:,ii) = computeR2(DATA,MODEL);
        
        % For specific conditions
        cond = unique(stimcond); nCond = length(cond);
        for jj = 1:nCond
            condInx = find(stimcond == cond(jj));
            DATA = flatten(data(:,condInx,ii));
            MODEL = flatten(pred(:,condInx,ii));
            R2cond(jj,ii) = computeR2(DATA,MODEL);
        end
        
        %fprintf('[%s] R2 for dataset %d = %0.2f \n', mfilename, ii, R2concat(:,ii))
    
        %% COMPUTE MODELBASED DERIVED PARAMETERS 
        
        if includeDerivedParams
            % Compute parameters and generate a prediction to a sustained stimulus:
            [derived_prm, pred_names, pred_sustained, pred_transient] = tde_computeDerivedParams(objFunction, params(:,ii));

            derivedPrm(:,ii) = derived_prm;
            derivedPredS(:,ii) = pred_sustained;
            derivedPredT(:,ii) = pred_transient;
        end
        
    end
    
    %% COLLECT RESULTS
    results(kk).model           = objFunction;
	results(kk).params          = params;
    results(kk).pred            = pred;
    results(kk).R2.stim         = R2stim;
    results(kk).R2.concat_all   = R2concat;
    results(kk).R2.concat_cond  = R2cond; 
    results(kk).derived.names   = pred_names;
    results(kk).derived.params  = derivedPrm;
    results(kk).derived.pred_s  = derivedPredS;
    results(kk).derived.pred_t  = derivedPredT;
    
end
fprintf('[%s] Done!\n', mfilename)
end

