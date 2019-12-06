function [results] = tde_evaluateModelFit(objFunction, params, data, pred, makePlots)


nDatasets   = size(data,3);
nStim       = size(data,2);

derivedPrm  = nan(2,nDatasets);
rSq         = nan(nStim,nDatasets);

%% EXTRACT SUMMARY METRICS 

for ii = 1:nDatasets % loop over channels or channel averages

    [derived_prm, pred_derived] = tde_computeDerivedParams(objFunction, params(:,ii));

    derivedPrm(1,ii) = derived_prm.t2pk;
    derivedPrm(2,ii) = derived_prm.r_asymp;
    derivedPred(:,ii)= pred_derived;
    
    rsq = nan(nStim,1);
    for jj = 1:nStim % loop over stimuli
        mdl = fitlm(pred(:,jj,ii), data(:,jj, ii));
        rsq(jj) = mdl.Rsquared.Ordinary;
    end
    disp(mean(rsq));
    rSq(:,ii) = rsq;
    % TO DO: compute RSQ for concatenated data as well (per condition and
    % all)
end

results.rSquare     = rSq;
results.fittedPrm   = params;
results.derivedPrm  = derivedPrm;
results.derivedPred = derivedPred;

return
%% MAKE PLOTS
if makePlots
    
%% plot fitted and derived Params

figure;hold on
derivedTitles = {'time2peak', 'R_asymp'};
fittedTitles = {'tau1', 'weight','tau2', 'n', 'sigma'};
nChans = height(channels);
subplot(2,5,1); plot(1:nChans,mean(results.rSquare), '.b', 'MarkerSize', 50, 'LineStyle', 'none')
set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
title('explained variance'); xlabel('visual area');  ylabel('R2'); set(gca, 'fontsize', 16);

for p1 = 1:2
    subplot(2,5,p1+1);
    plot(1:nChans,results.derivedPrm(p1,:), '.r', 'MarkerSize', 50, 'LineStyle', 'none')
    set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
    title(derivedTitles{p1}); xlabel('visual area');  ylabel('parameter value'); set(gca, 'fontsize', 16);
end
for p1 = 1:5
    subplot(2,5,p1+5);
    plot(1:size(data2fit,3),results.fittedPrm(p1,:), '.k', 'MarkerSize', 50, 'LineStyle', 'none')
    set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
    title(fittedTitles{p1}); xlabel('visual area'); ylabel('parameter value');set(gca, 'fontsize', 16);
end
set(gcf, 'Position', [400 200 2000 1200]);


%% plot timecourses and predictions

% visualization 2
conditionsOfInterest = {'CRF', 'ONEPULSE', 'TWOPULSE'};
nCond = length(conditionsOfInterest);
for kk = 1:size(smallData,3)
    figure('Name', channels.name{kk});
    d = smallData(:,:,kk);
    maxresp = max(d(:));
    p1 = pred1(:,:,kk);
    p2 = pred2(:,:,kk);
    for ii = 1:length(conditionsOfInterest)
        inx = contains(stimnames, conditionsOfInterest{ii});
        subplot(nCond,1,ii); hold on
        h = plot(flatten(stim_ts(:,inx))*maxresp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        plot(flatten(d(:,inx)), 'k', 'LineWidth', 2); 
        plot(flatten(p1(:,inx)), 'b', 'LineWidth', 2);
        plot(flatten(p2(:,inx)), 'r', 'LineWidth', 2);
        legend('data', 'TTCSTIG19','DN'); 
        title(sprintf('%s r2 DN = %0.2f , r2 TTCSTIG19 = %0.2f',conditionsOfInterest{ii}, mean(results2.rSquare(inx,kk)), mean(results1.rSquare(inx,kk))));
        axis tight   
        set(gca, 'XTick',1:length(t):length(find(inx))*length(t), 'XTickLabel', []);
    end
    set(gcf, 'Position', [400 200 1800 1200]);
end

%% plot prediction for sustained stimulus

figure;hold on
plot(results1.derivedPred, 'LineWidth', 2);
plot(results2.derivedPred, 'LineWidth', 2);
set(gca, 'Xlim', [0 1000]);
leg1 = sprintf('TTCSTIG t2p = %0.2f rasymp = %0.2f', results1.derivedPrm(1), results1.derivedPrm(2));
leg2 = sprintf('DN t2p = %0.2f rasymp = %0.2f', results2.derivedPrm(1), results2.derivedPrm(2));
legend({leg1,leg2}); 
set(gca, 'FontSize', 12);
title('Model predictions for sustained stimulus');

end

