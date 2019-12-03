function [results] = tde_evaluateModelFit(objFunction, params, data, pred)


nDatasets   = size(data,3);
nStim       = size(data,2);

derivedPrm  = nan(2,nDatasets);
rSq         = nan(nStim,nDatasets);

%% EXTRACT SUMMARY METRICS 

for ii = 1:nDatasets % loop over channels or channel averages

    %derived_prm = tde_computeDerivedParams(objFunction, params(:,ii));

    %derivedPrm(1,ii) = derived_prm.t2pk;
    %derivedPrm(2,ii) = derived_prm.r_asymp;
    
    r = nan(nStim,1);
    for jj = 1:nStim % loop over stimuli
        mdl = fitlm(pred(:,jj,ii), data(:,jj, ii));
        r(jj) = mdl.Rsquared.Ordinary;
    end
    disp(mean(r));
    rSq(:,ii) = r;
end

results.rSquare     = rSq;
results.derivedPrm  = derivedPrm;
results.fittedPrm   = params;

