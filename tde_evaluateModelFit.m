function [results] = tde_evaluateModelFit(data, objFunction, params, pred)

% Inputs: cell arrays with objFunction ,params, data, pred
if ~iscell(objFunction), objFunction = {objFunction}; end
if ~iscell(params), objFunction = {params}; end
if ~iscell(pred), objFunction = {pred}; end
% if ~exist(makePlots, 'var') || isempty(makePlots), makePlots = 0; end

nModels     = size(objFunction,2);
nStim       = size(data,2);
nDatasets   = size(data,3);

rSq_bystim  = nan(nStim,nDatasets);
rSq_concat  = nan(1,nDatasets);
derivedPrm  = nan(2,nDatasets);

%% COMPUTE SUMMARY METRICS
results = [];

% Loop over models
for kk = 1:nModels
    
    
    % Loop over channels or channel averages
    for ii = 1:nDatasets
        
        %% COMPUTE R-SQUARE
        
        % For each individual stimulus 
        rsq = nan(nStim,1);
        for jj = 1:nStim % loop over stimuli
            mdl = fitlm(pred{kk}(:,jj,ii), data(:,jj, ii));
            rsq(jj) = mdl.Rsquared.Ordinary;
        end
        rSq_bystim(:,ii) = rsq;
        
        % For all stimuli concatenated
        mdl = fitlm(flatten(pred{kk}(:,:,ii)), flatten(data(:,:,ii)));
        rSq_concat = mdl.Rsquared.Ordinary;
    end

    %% COMPUTE DERIVED PARAMETERS 
    for ii = 1:nDatasets % loop over channels or channel averages
        
        % Generate a prediction to a sustained stimulus:
        [derived_prm, pred_derived] = tde_computeDerivedParams(objFunction{kk}, params{kk}(:,ii));

        derivedPrm(1,ii) = derived_prm.t2pk;
        derivedPrm(2,ii) = derived_prm.r_asymp;
        derivedPred(:,ii)= pred_derived;
    end

    %% COLLECT RESULTS
    results(kk).model       = objFunction{kk};
	results(kk).params      = params{kk};
    results(kk).pred        = pred{kk};
    results(kk).rSquareStim = rSq_bystim;
    results(kk).rSquareConc = rSq_concat;
    results(kk).derivedPrm  = derivedPrm;
    results(kk).derivedPred = derivedPred;
end

end

