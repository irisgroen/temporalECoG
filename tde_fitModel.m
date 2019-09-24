function [out] = tde_fitModel(data, stims_ts, functionhandle, opts, savePlots)
% Description

	

fittedPrm = [];
derivedPrm = [];

    %% scale each electrode to the max (make optional?)
    normdata = [];
    for k = 1 : size(data,1)
        tmp = data(k, :, :);
        maxRsp(k) = max(tmp(:));
        normdata(k, :, :) = data(k, :, :)./maxRsp(k);
    end
    
    %% average elecs within area (make optional)  
    INX = [];
    INX{1} = contains(elec_info.wang, 'V1') | contains(elec_info.benson, 'V1');
    INX{2} = contains(elec_info.wang, 'V2') | contains(elec_info.benson, 'V2');
    INX{3} = contains(elec_info.wang, 'V3') & ~contains(elec_info.wang, {'V3a', 'V3b'}) | contains(elec_info.benson, 'V3') & ~contains(elec_info.benson, {'V3a', 'V3b'});
    INX{4} = contains(elec_info.wang, {'V3a', 'V3b', 'IPS', 'hV4', 'LO', 'TO'}) | contains(elec_info.benson, {'V3a', 'V3b', 'hV4', 'LO', 'TO'}) ;

     switch whichArea
        case 'V1'
            data2fit = normdata(INX{1},:,:); % V1 ELECS
        case 'V2'
            data2fit = normdata(INX{2},:,:); % V2 ELECS
        case 'V3'
            data2fit = normdata(INX{3},:,:); % V3 ELECS
        case 'Vhigher'
            data2fit = normdata(INX{4},:,:); % higher ELECS
     end
    
     % or scale only here?
    mdata    = squeeze(mean(data2fit,1));
    maxmdata = max(mdata(:));
    mdata    = mdata ./maxmdata;
    
    %% FIT THE DN model

    seed = [0.03, 0.07, 1.5, 0.15, 0.06, 1];
    %seed = [0.1, 0.1, 3, 0.1, 0.06, 1];
    lb   = [0, 0, 0, 0, 0, 0];
    ub   = [1, 1, 10, 1, 1, 1];

    prm = [];
    prm = fminsearchbnd(@(x) dn2_fineFitCtrstDur(x, mdata, t, stim), seed, lb, ub);
    
    %% GENERATE MODEL PREDICTIONS

    prm_tofit = [prm(1), 0, prm(2 : end)];

    pred = dn_DNmodel(prm_tofit, stim, t);
    pred = pred./max(pred(:));
    
    %% EXTRACT SUMMARY METRICS 
   
	derived_prm = dn_computeDerivedParams(prm, 'uniphasic');
    
    derivedPrm(aa,1) = derived_prm.t2pk;
    derivedPrm(aa,2) = derived_prm.r_asymp;
    fittedPrm(aa,:) = prm;
    
    
     %% CALCULATE R2
    rvals = [];
    for k = 1:17
        rvals(k) = corr(pred(k,timeInd)',mdata(k,timeInd)');  
    end
    mean(rvals.^2)
    
    %% VARIOUS PLOTS from old code
    
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