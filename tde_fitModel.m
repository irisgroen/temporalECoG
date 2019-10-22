function [results, pred] = tde_fitModel(objFunction, data, stim, opts)
% Description
% To call fitting function, we need:
%   1. the objective function (model form and type of error)
%   2. data
%   3. stimuli
%   4. starting values and bounds for parameters (there should be defaults
%               for each type of model)
%   5. sample rate of the data
 
%% FIT THE temporal model

% model start point and bounds
x0 = opts.x0;
lb = opts.lb;
ub = opts.ub;

if isfield(opts, 'pub')
    plb = opts.plb;
    pub = opts.pub;
end

% sample rate (Hz)
srate = opts.srate;

% initialize
nParams     = length(x0);
nDatasets   = size(data,3);
nStim       = size(stim,2);

fittedPrm   = nan(nParams,nDatasets);
derivedPrm  = nan(2,nDatasets);
rSq         = nan(nStim,nDatasets);

options = optimset('Display','iter');

% options.MaxFunEvals = 10000;

for ii = 1:nDatasets % loop over channels or channel averages

    fprintf('[%s] Fitting model for dataset %d \n',mfilename, ii);
    
    data2fit = data(:,:,ii);
    
%     prm = fminsearchbnd(@(x) objFunction(x, data2fit, stim, srate), x0, lb, ub, options);
%     
%     % test bads by putting in a good set of initial parameters
%     x0 = prm;
    
    prm =bads(@(x) objFunction(x, data2fit, stim, srate),  x0, lb, ub, plb, pub, [], options);
                
    
    %% GENERATE MODEL PREDICTIONS

    [~, pred] = objFunction(prm, [], stim, srate);

    % pred = pred./max(pred(:));

    %% EXTRACT SUMMARY METRICS 

    derived_prm = tde_computeDerivedParams(objFunction, prm);

    derivedPrm(1,ii) = derived_prm.t2pk;
    derivedPrm(2,ii) = derived_prm.r_asymp;
    fittedPrm(:,ii) = prm;

    %% CALCULATE R2
    r = nan(nStim,1);
    for jj = 1:nStim % loop over stimuli
        mdl = fitlm(pred(:,jj), data2fit(:,jj));
        r(jj) = mdl.Rsquared.Ordinary;
    end
    disp(mean(r));
    rSq(:,ii) = r;
end

results.derivedPrm  = derivedPrm;
results.fittedPrm   = fittedPrm;
results.rSquare     = rSq;

return    
    
    
    %% VARIOUS PLOTS from old code (need to be updated or moved to separate function)
    
    %% VISUALIZE ELECTRODE RESPONSES
    figure (1); clf;
    
    v = find(contains(AreaNames, whichArea));
    inx = find(INX{v});
    %colors = jet(length(inx));
    nSubplot = ceil(sqrt(length(inx))); names = [];
    for ii  = 1:length(find(inx))
        %subplot(nSubplot,nSubplot,ii)
        patientIdx = find(contains(channels.subject_name, elec_info.patientname{inx(ii)}));
        ecog_plotSingleTimeCourse(t,squeeze(mean(data(inx(ii),:,:),2)), [], colors(patientIdx,:));
        %names{ii} = sprintf('%s %s %s W:%s B:%s', elec_info.patientname{inx(ii)}, elec_info.sessionname{inx(ii)}, elec_info.elecname{inx(ii)}, elec_info.wang{inx(ii)},elec_info.benson{inx(ii)});
    end
    legend(names)
    title(sprintf('%s (n = %d)', AreaNames{v}, length(find(INX{v}))));
    set(1, 'Position', [800 800 900 900]);
    
    %% VISUALIZE MODEL FIT
    figure (2), clf
    for k = 1 : 17
       subplot(17, 1, k)
       plot(t, mdata(k, :), 'r-', 'linewidth', 2), hold on, axis tight, ylim([-0.1, 1]), box off
    end
    
    figure (2)
    for k = 1 : 17
       subplot(17, 1, k)
       plot(t, pred(k, :), 'k-', 'linewidth', 2)
       set(gca, 'ytick', [0, 1], 'xtick', [0, 0.5, 1]), xlim([-0.2, 1])
       if k == 17, xlabel('time (s)'), end
    end
    set(2, 'Position', [100 100 500 1300]);

    %% ALTERNATIVE WAY TO VISUALIZE DATA

    figure (3), clf
    subplot(2, 3, 1), set(gca, 'colororder', cool(5)), hold on
    plot(t, mdata(1 : 5, :), '-', 'linewidth', 2), 
    subplot(2, 3, 4), set(gca, 'colororder', cool(5)), hold on
    plot(t, pred(1 : 5, :), '-', 'linewidth', 2), 

    subplot(2, 3, 2), set(gca, 'colororder', cool(6)), hold on
    plot(t, mdata(6 : 11, :), '-', 'linewidth', 2), 

    subplot(2, 3, 5), set(gca, 'colororder', cool(6)), hold on
    plot(t, pred(6 : 11, :), '-', 'linewidth', 2), 

    subplot(2, 3, 3), set(gca, 'colororder', cool(6)), hold on
    plot(t, mdata(12 : 17, :), '-', 'linewidth', 2), 
    subplot(2, 3, 6), set(gca, 'colororder', cool(6)), hold on
    plot(t, pred(12 : 17, :), '-', 'linewidth', 2),

    for k = 1 : 6
        subplot(2, 3, k)
        axis tight, xlim([-0.2, 1]), set(gca, 'ytick', [0, 1], 'xtick', [0, 0.5, 1], 'fontsize', 16), ylim([0, 1])
        xlabel('time (s)')
    end

    %% plot SUM of prediction

    timeInd = t>0 &t<1;
    sumData = sum(mdata(:,timeInd),2);
    sumPred = sum(pred(:,timeInd),2);
    
    figure (4), clf;

    % PREDICTIONS
    
    % compute linpred for contrast  
	Y = sumPred(1:5);
    X = contrasts';
    beta = regress(Y,X);   
    linpred = contrasts*beta;
    
    subplot(2,3,4); hold on;
    plot(contrasts, sumPred(1:5), 'Marker', '.', 'MarkerSize', 50,'LineStyle', '-', 'LineWidth', 3, 'Color', 'r');
    plot(contrasts, linpred,'LineStyle', '-','Color','g', 'LineWidth',3);
    xlabel('contrast (%)'); ylabel('sum (0-1s)'); title('model prediction');
    axis tight, xlim([0, 1.2]), set(gca, 'fontsize', 16), ylim([0, max([sumPred' linpred])+0.1*max([sumPred' linpred])])
    legend({'DN model', 'linear prediction'}, 'Location', 'NorthWest');
    set(gca, 'XTick', contrasts, 'XTickLabel', round(contrasts,3), 'XTickLabelRotation', 45);box off
    
    % compute linpred for duration  
	Y = sumPred(6:11);
    X = durs';
    beta = regress(Y,X);   
	linpred = durs*beta;

    subplot(2,3,5); hold on;
    plot(durs, sumPred(6:11), 'Marker', '.', 'MarkerSize', 50,'LineStyle', '-', 'LineWidth', 3, 'Color', 'r');
    plot(durs, linpred,'LineStyle', '-','Color','g', 'LineWidth',3);
    xlabel('stimulus duration (s)'); ylabel('sum (0-1s)'); title('model prediction');
    axis tight, xlim([0, 0.6]), set(gca, 'fontsize', 16), ylim([0, max([sumPred' linpred])+0.1*max([sumPred' linpred])])
    set(gca, 'XTick', durs, 'XTickLabel', round(durs,3), 'XTickLabelRotation', 45);box off

    subplot(2,3,6), hold on;
    plot(stimISIs, sumPred([10 12:17]), 'Marker', '.', 'MarkerSize', 50,'LineStyle', '-','LineWidth', 3, 'Color', 'r');
    plot(stimISIs, ones(length(stimISIs))*durs(5)*beta,'LineStyle', '-','Color','g', 'LineWidth',3);
    xlabel('stimulus ISI (s)'); ylabel('sum (0-1s)');title('model prediction');
    axis tight, xlim([0, 0.6]), set(gca, 'fontsize', 16), ylim([0, max([sumPred' linpred])+0.1*max([sumPred' linpred])])
    set(gca, 'XTick', stimISIs, 'XTickLabel', round(stimISIs,3), 'XTickLabelRotation', 45);box off
    
    % DATA
    
    % compute linpred for contrast  
	Y = sumData(1:5);
    X = contrasts';
    beta = regress(Y,X);   
    linpred = contrasts*beta;   
    
    subplot(2,3,1); hold on;
    plot(contrasts, sumData(1:5), 'Marker', '.', 'MarkerSize', 50,'LineStyle', '-', 'LineWidth', 3, 'Color', 'k');
    plot(contrasts, linpred,'LineStyle', '-','Color','g', 'LineWidth',3);
    xlabel('contrast (%)'); ylabel('sum (0-1s)'); title('ecog data');
    axis tight, xlim([0, 1.2]), set(gca, 'fontsize', 16), ylim([0, max([sumData' linpred])+0.1*max([sumData' linpred])])
    legend({'ecog data', 'linear prediction'}, 'Location', 'NorthWest');
    set(gca, 'XTick', contrasts, 'XTickLabel', round(contrasts,3)*100, 'XTickLabelRotation', 45);box off
    
    % compute linpred for duration  
	Y = sumData(6:11);
    X = durs';
    beta = regress(Y,X);   
	linpred = durs*beta;
    
    subplot(2,3,2); hold on;
    plot(durs, sumData(6:11), 'Marker', '.', 'MarkerSize', 50,'LineStyle', '-', 'LineWidth', 3, 'Color', 'k');
    plot(durs, linpred,'LineStyle', '-','Color','g', 'LineWidth',3);
    xlabel('stimulus duration (s)'); ylabel('sum (0-1s)'); title('ecog data');
    axis tight, xlim([0, 0.6]), set(gca, 'fontsize', 16), ylim([0, max([sumData' linpred])+0.1*max([sumData' linpred])])
    set(gca, 'XTick', durs, 'XTickLabel', round(durs,3), 'XTickLabelRotation', 45);box off

    subplot(2,3,3), hold on;
    plot(stimISIs, sumData([10 12:17]), 'Marker', '.', 'MarkerSize', 50,'LineStyle', '-','LineWidth', 3, 'Color', 'k');
    plot(stimISIs, ones(length(stimISIs))*durs(5)*beta,'LineStyle', '-','Color','g', 'LineWidth',3);
    xlabel('stimulus ISI (s)'); ylabel('sum (0-1s)');title('ecog data');
    axis tight, xlim([0, 0.6]), set(gca, 'fontsize', 16), ylim([0, max([sumData' linpred])+0.1*max([sumData' linpred])])
    set(gca, 'XTick', stimISIs, 'XTickLabel', round(stimISIs,3), 'XTickLabelRotation', 45);box off

    set(4, 'Position', [100 200 1200 900]);

   
    %% plot fitted parameters
    figure (7) ;hold on
    derivedTitles = {'time2peak', 'R_asymp'};
    fittedTitles = {'tau1', 'tau2', 'n', 'sigma'};

    for p = 1:2
        subplot(2,4,p);
        plot(1:4,derivedPrm(:,p), '.r', 'MarkerSize', 50, 'LineStyle', 'none')
        set(gca, 'Xlim', [0 5], 'XTick', 1:4, 'XTickLabel', AreaNames);
        set(gca, 'fontsize', 16);
        title(derivedTitles{p}); xlabel('visual area');  ylabel('parameter value');
    end


    for p = 1:4
        subplot(2,4,p+4);
        plot(1:4,fittedPrm(:,p), '.k', 'MarkerSize', 50, 'LineStyle', 'none')
        set(gca, 'Xlim', [0 5], 'XTick', 1:4, 'XTickLabel', AreaNames);
        set(gca, 'fontsize', 16);
        title(fittedTitles{p}); xlabel('visual area'); ylabel('parameter value');
    end


    
    
end