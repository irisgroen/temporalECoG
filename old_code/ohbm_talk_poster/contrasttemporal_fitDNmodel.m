% ADD PATHS
addpath(genpath('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal'));
addpath(genpath('/Volumes/server/Projects/Temporal_integration/DN_2018_data/code'))

% LOAD ECOG DATA
fileName = 'data_contrasttemporal.mat';
load(fullfile('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess', fileName));

% LOAD FMRI DATA
load('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess/fMRIdataZhouPCB')
load('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess/preproc_BAIRfmri.mat')

% LOAD STIMULUS DATA 
stimData = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/stimuli/sub-som723_ses-nyuecog02_task-temporalpattern_run-01.mat');
stimDurations = unique(stimData.stimulus.duration);
stimISIs = unique(stimData.stimulus.ISI);
    
stimData = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/stimuli/sub-som723_ses-nyuecog02_task-spatialpattern_run-01.mat');
CRFstimuli = stimData.stimulus.im_cell(contains(stimData.stimulus.categories, 'CRF'));
CRFrms = []; CRFp2p = [];
for ii = 1:length(CRFstimuli)
    CRFrms(ii) = rms(double(CRFstimuli{ii}(:))/255-0.5);
    CRFp2p(ii) = peak2peak(double(CRFstimuli{ii}(:))/255-0.5);
end

% PREDEFINED BY JING
% contrasts = [0.0405, 0.0902, 0.2105, 0.3203, 1.0000];
% durs      = [0.016667, 0.033333, 0.066667, 0.13333, 0.26667, 0.53333];

% PREDEFINED IN STIM GENERATION
contrastLevels = [0.0625 0.125 0.25 0.5 1]; 

contrasts = contrastLevels;
durs = stimDurations;

% DEFINE AREAS
AreaNames = {'V1', 'V2', 'V3', 'Vhigher'};
patientList = unique(elec_info.patientname);
colors = parula(length(patientList));

% PLOT SPECS
saveFigure = 0;
figureFormat = 'EPS';  % PNG or EPS
figureFolder = 'nyu_and_umcu_thresh1.5'; % 

%% shift chaam elecs

% exclude chaam elecs
% tokeep = contains(elec_info.patientname, 'umcu');
% data = data(tokeep,:,:);
% elec_info = elec_info(tokeep,:);
% 
% shiftInSamples = 32;
% 
% ind = find(contains(elec_info.patientname, 'chaam'));
% toshift = data(ind,:,:);
% newdata = zeros(size(toshift));
% newdata(:,:,1:size(toshift,3)-shiftInSamples) = toshift(:,:,shiftInSamples+1:end);
% newdata(:,:,size(newdata,3)+shiftInSamples:end) = 0;
% 
% data(ind,:,:) = newdata;

%%
allPrm = [];
derivedPrm = [];

for aa = 1:1%length(AreaNames)
    
    whichArea = AreaNames{aa};
    
    %% Compute average by area

    INX = [];
    INX{1} = contains(elec_info.wang, 'V1') | contains(elec_info.benson, 'V1');
    INX{2} = contains(elec_info.wang, 'V2') | contains(elec_info.benson, 'V2');
    INX{3} = contains(elec_info.wang, 'V3') & ~contains(elec_info.wang, {'V3a', 'V3b'}) | contains(elec_info.benson, 'V3') & ~contains(elec_info.benson, {'V3a', 'V3b'});
    INX{4} = contains(elec_info.wang, {'V3a', 'V3b', 'IPS', 'hV4', 'LO', 'TO'}) | contains(elec_info.benson, {'V3a', 'V3b', 'hV4', 'LO', 'TO'}) ;

    %% PLOT average response for each electrode
    figure (1); clf;
    
    v = find(contains(AreaNames, whichArea));
    inx = find(INX{v});
    %colors = jet(length(inx));
    nSubplot = ceil(sqrt(length(inx))); names = [];
    for ii  = 1:length(find(inx))
        %subplot(nSubplot,nSubplot,ii)
        patientIdx = find(contains(patientList, elec_info.patientname{inx(ii)}));
        ecog_plotSingleTimeCourse(t,squeeze(mean(data(inx(ii),:,:),2)), [], colors(patientIdx,:));
        names{ii} = sprintf('%s %s %s W:%s B:%s', elec_info.patientname{inx(ii)}, elec_info.sessionname{inx(ii)}, elec_info.elecname{inx(ii)}, elec_info.wang{inx(ii)},elec_info.benson{inx(ii)});
    end
    legend(names)
    title(sprintf('%s (n = %d)', AreaNames{v}, length(find(INX{v}))));
    set(1, 'Position', [800 800 900 900]);
    
    %% COMPUTE THE MAXIMUM RESPONSE IN EACH SUBJECT'S DATA SET
    normdata = [];
    for k = 1 : size(data,1)
        tmp = data(k, :, :);
        maxRsp(k) = max(tmp(:));
        normdata(k, :, :) = data(k, :, :)./maxRsp(k);
    end

    normdata = data;
    
    %% AVERAGE THE RESPONSE ACROSS SUBJECTS AND ELECTRODES

    switch whichArea
        case 'V1'
            data_to_fit = normdata(INX{1},:,:); % V1 ELECS
        case 'V2'
            data_to_fit = normdata(INX{2},:,:); % V2 ELECS
        case 'V3'
            data_to_fit = normdata(INX{3},:,:); % V3 ELECS
        case 'Vhigher'
            data_to_fit = normdata(INX{4},:,:); % higher ELECS
    end

    mdata    = squeeze(mean(data_to_fit,1));
    maxmdata = max(mdata(:));
    mdata    = mdata ./maxmdata;

    figure (2), clf
    for k = 1 : 17
       subplot(17, 1, k)
       plot(t, mdata(k, :), 'r-', 'linewidth', 2), hold on, axis tight, ylim([-0.1, 1]), box off
    end

    %% MAKE STIMULUS

    nStim = size(data_to_fit, 2);
    stim  = zeros(nStim, length(t));

    % CONTRAST STIMULI --------------------------------------------------------
    stim(1 : 5, t > 0 & t<=0.5) = 1;
    for k = 1 : 5, stim(k, :) = stim(k, :) .* contrasts(k); end

    % INCREASING DURATIONS ----------------------------------------------------
    for k = 1 : 6, stim(k + 5, (t>0) & (t <= durs(k))) = 1; end

    % INCREASING ISI ----------------------------------------------------------
    stim(12 : nStim, t > 0 & t <= durs(4)) = 1;
    for k = 1 : 6
        t_start = durs(4) + durs(k);
        t_end   = durs(4) * 2 + durs(k);
        stim(11 + k, t > t_start & t <= t_end) = 1;
    end

    % VISUALIZE THE STIMULI ---------------------------------------------------
    figure (2)
    for k = 1 : 17
       subplot(17, 1, k)
       plot(t, stim(k, :), 'k-')
    end

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
    
    %% VISUALIZE MODEL FIT

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

    %% CALCULATE R2
    rvals = [];
    for k = 1:17
        rvals(k) = corr(pred(k,timeInd)',mdata(k,timeInd)');  
    end
    mean(rvals.^2)

    %% PLOT SUM vs FMRI ZHOU (code copied from /Volumes/server/Projects/Temporal_integration/DN_2018_code_data/code/dn_fitDNECoG2fMRI.m'
    
    roiNm = {'V1', 'V2', 'V3'};

    linrsp = [1, 2, 4, 8, 16, 32, 16, 16, 16, 16, 16, 16, 0];
    x1 = [0, 1, 2, 4, 8, 16, 32];
    
    mriOrder{1} = [13,1:6]; % onepulse
    mriOrder{2} = [5 7:12]; % twopulse
    ecogOrder{1} = [6:11]; % onepulse
    ecogOrder{2} = [10 12:17]; % twopulse
    
    switch whichArea
        case 'V1'
            iroi = 1;
            %ecogScaleFac = 0.0035;
        case 'V2'
            iroi = 2;
            %ecogScaleFac = 0.0035;
        case 'V3'
            iroi = 3;
            %ecogScaleFac = 0.0035;
        otherwise
            iroi = [];
    end

    if ~isempty(iroi)
        figure (5); clf    
        
        % compute scale factor    
        Y = mfmri_2fit(iroi, [mriOrder{1} mriOrder{2}]);
        X = [0 sumPred([ecogOrder{1} ecogOrder{2}])'];
        ecogScaleFacPred = regress(Y',X');
        X = [0 sumData([ecogOrder{1} ecogOrder{2}])'];
        ecogScaleFacData = regress(Y',X');
        
        for jj = 1:2
           
            subplot(1,2,jj);
            %subplot(2, 3, iroi + 3*(jj-1))
            %figure;hold on;
            
            plot(x1, mfmri_2fit(iroi, mriOrder{jj}), 'ko', 'markersize', 9,  'markerfacecolor', 'k'), hold on
            
            for k = 1 : length(mriOrder{jj})
                low = sfmri_2fit(iroi, 1, mriOrder{jj}(k));
                high = sfmri_2fit(iroi, 2, mriOrder{jj}(k));
                p = plot(x1(k)*[1 1], [low, high], '-', 'linewidth', 3, 'color', 0.8 * ones(1, 3));
                p.Annotation.LegendInformation.IconDisplayStyle = 'off';
            end
            %xlim([0.5, length(mriOrder) + 0.5]),
            %ylim([-.05, 0.6]), xlim([-2, 34]), box off, set(gca, 'xtick', [0, 4, 16, 32], 'xticklabel', [0, 67, 267, 533])
            
            title(roiNm{iroi}), set(gca, 'fontsize', 16)
            
            plot(x1, linscale(iroi).*linrsp(mriOrder{jj}), 'g-', 'linewidth', 3)
            plot(x1, ecogPrd1(iroi, mriOrder{jj}), 'r-', 'linewidth', 3)
            if jj == 1
                plot(x1, [0; sumPred(ecogOrder{jj})*ecogScaleFacPred], 'c-', 'linewidth', 3)
                plot(x1, [0; sumData(ecogOrder{jj})*ecogScaleFacData], 'm-', 'linewidth', 3)
                legend({'BOLD', 'linear pred', 'ECOG fit PCB', 'ECOG fit', 'ECOG data'}, 'Location', 'NorthWest');
                set(gca, 'XTick', x1, 'XTickLabel', [0 round(stimDurations,3)],  'XTickLabelRotation', 45);
            else
                plot(x1, sumPred(ecogOrder{jj})*ecogScaleFacPred, 'c-', 'linewidth', 3)
                plot(x1, sumData(ecogOrder{jj})*ecogScaleFacData, 'm-', 'linewidth', 3)
                set(gca, 'XTick', x1, 'XTickLabel', round(stimISIs,3), 'XTickLabelRotation', 45);
            end
            set(gca, 'xaxislocation', 'origin'), xlabel('time (ms)'), ylabel('% BOLD')
        end
        box off
    end

    set(5, 'Position', [150 100 1500 700]);
    
    %% PLOT SUM vs FMRI BAIR
    
    %linrsp = [1, 2, 4, 8, 16, 32, 16, 16, 16, 16, 16, 16, 0];
    %x1 = [0, 1, 2, 4, 8, 16, 32];
    
    condOrder{1} = 1:5; % CRF
    condOrder{2} = 6:11; %ONEPULSE
    condOrder{3} = [10 12:17]; %TWOPULSE
    
    figure (6); clf 
    
    if aa > 3
        mB = squeeze(mean(mean(B(:,aa:end,:),2),1)); %Vhigher
        mB_se = squeeze(std(mean(B(:,aa:end,:),2),0,1))/sqrt(size(B,1)); %Vhigher
    else
        mB = squeeze(mean(B(:,aa,:),1)); %V1,2,3
        mB_se = squeeze(std(B(:,aa,:),0,1))/sqrt(size(B,1));
    end
    
    % compute scale factor  
    Y = mB;
    X = sumPred;
    ecogScaleFacPred = regress(Y,X);
    X = sumData;
    ecogScaleFacData = regress(Y,X);
   
    for jj = 1:3
        
        if jj == 1
            x1 = contrasts*100; 
            predScaleFac = regress(Y(condOrder{jj}),x1');
            linrsp = x1*predScaleFac; 
            plot_title = [whichArea ' contrast']; 
            x_unit = 'contrast (%)';
            xtickloc = contrasts*100; xtickvals = round(contrasts,3)*100;
            subplot(1,3,3);
        elseif jj == 2
            x1 = durs;
            predScaleFac = regress(Y([condOrder{jj}]),x1');
            linrsp = x1*predScaleFac;
            plot_title = [whichArea ' duration']; 
            x_unit = 'time (ms)';
            xtickloc = durs; xtickvals = round(durs,3);
            subplot(1,3,1);
        elseif jj == 3
            x1 = stimISIs;
            linrsp = ones(length(x1),1) * durs(5)*predScaleFac;
            plot_title = [whichArea ' ISI']; 
            x_unit = 'time (ms)';
            xtickloc = stimISIs; xtickvals = round(stimISIs,3);
            subplot(1,3,2);
        end
    
        plot(x1, mB(condOrder{jj}), 'ko', 'markersize', 20,  'markerfacecolor', 'k'), hold on

        for k = 1 : length(condOrder{jj})
            low = mB(condOrder{jj}(k)) - mB_se(condOrder{jj}(k));%sfmri_2fit(iroi, 1, mriOrder{jj}(k));
            high = mB(condOrder{jj}(k)) + mB_se(condOrder{jj}(k));
            p = plot(x1(k)*[1 1], [low, high], '-', 'linewidth', 3, 'color', 0.5 * ones(1, 3));
            p.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
               
        % plot linear prediction
        plot(x1,linrsp, 'k:', 'linewidth', 3)
        plot(x1, sumPred(condOrder{jj})*ecogScaleFacPred, 'g-', 'Marker', '.', 'MarkerSize', 50,'linewidth', 3)
        plot(x1, sumData(condOrder{jj})*ecogScaleFacData, 'm-', 'Marker', '.', 'MarkerSize', 50,'linewidth', 3)
            
        if jj == 1
            %legend({'BOLD', 'linear pred', 'ECOG fit', 'ECOG data'});
            %legend({'BOLD', 'linear pred', 'ECOG DN model fit'});
        end
        set(gca, 'YLim', [0 1.2], 'XLim', [min([0 x1])-0.1*min([0 x1]) max(x1)+0.1*max(x1)]);
        set(gca, 'xaxislocation', 'origin'), xlabel(x_unit), ylabel('% BOLD');
        set(gca, 'XTick', xtickloc, 'XTickLabel', xtickvals, 'XTickLabelRotation', 45);box off
        set(gca, 'fontsize', 16);box off;
        %title(plot_title);
    end
    
    %set(6, 'Position', [250 250 1500 700]);
    set(6, 'Position', [200 200 2000 600]);
    %% SAVE FIGURES

    if saveFigure

        fg1Nm = sprintf('broadbandByElectrodes_%s', whichArea);
        fg2Nm = sprintf('visualizeModelFit1_%s', whichArea);
        fg3Nm = sprintf('visualizeModelFit2_%s', whichArea);
        fg4Nm = sprintf('sumModelFits_%s', whichArea);
        fg5Nm = sprintf('comparisonWithfMRI_Zhou_%s', whichArea);
        fg6Nm = sprintf('comparisonWithfMRI_BAIR_%s', whichArea);

        saveLoc = fullfile(dn_ctrst_RootPath, 'analysisFigures',figureFolder);
        switch figureFormat
            case 'PNG'
                fmt = [1 300];
            case 'EPS'
                fmt = 0;
        end
        
        printnice(1, fmt, saveLoc, fg1Nm);
        printnice(2, fmt, saveLoc, fg2Nm);
        printnice(3, fmt, saveLoc, fg3Nm);
        printnice(4, fmt, saveLoc, fg4Nm);

        if ~isempty(iroi), printnice(5, fmt, saveLoc, fg5Nm); end
        
        printnice(6, fmt, saveLoc, fg6Nm);

    end
end

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


set(7, 'Position', [500 500 1500 800]);
 
if saveFigure
    fg7Nm = sprintf('summaryMetricsAndFittedModelParametersByArea', whichArea);
	printnice(7, fmt, saveLoc, fg7Nm);
end





