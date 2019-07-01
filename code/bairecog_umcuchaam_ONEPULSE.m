saveFigure = 0;

% Dataset specs
projectName = 'visual';
sub_label   = 'chaam'; 
ses_label   = 'UMCUECOGday03';

dataDir = '/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess';
dataName = fullfile(dataDir, sprintf('umcu%s_preproc_selectelecs',sub_label));

load(dataName);

% Check if we have the ECoG_utils repository on the path
if ~exist('ecog_plotTimecourses.m')
    tbUse ECoG_utils;
end

addpath(genpath('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal'));
addpath(genpath('/Volumes/server/Projects/Temporal_integration/DN_2018_code_data/code/external/'));

% EXTRACT STIMULUS DURATIONS
stimData = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/stimuli/sub-umcuchaam_ses-UMCUECOGday03_task-bairtemporalpattern_run-2_acq-clinical_events.mat');
stimDurations = unique(stimData.stimulus.duration);

eeToAverage = {[1:6],[6:8],[1:8]}; 
eeAvTitles = {'all V1 electrodes', 'all V2/V3 electrodes', 'all V1/V2/V3 electrodes'};

%% ONEPULSE %%
close all;

% PLOT time course for each condition
specs = [];
specs.dataTypes          = {'broadband'};
specs.smoothingLevelInMs = [];

specs.plot.colorMap      = 'parula';
specs.plot.nSubPlots     = [2 4];
specs.plot.addEccToTitle = 'yes';
specs.plot.showMax       = 'yes';

whichTrials = {'ONEPULSE-1', 'ONEPULSE-2','ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
[out] = ecog_plotTimecourses(trials, trials.channels.name, whichTrials, specs);

condName = whichTrials{1}(1:8);

% COMPUTE summary metrics (peak amplitude/sum/latency)
elecNames = trials.channels.name;
timeInx = trials.time>0.1 & trials.time<0.3;
elData = []; maxResp = []; maxRespSE = []; maxTmp = []; 

for ee = 1:length(elecNames)
    elData(:,:,ee) = out.broadband.(elecNames{ee}).mn;
    for ii = 1:length(whichTrials)
        [pks,locs] = findpeaks(elData(ii,timeInx,ee), trials.time(timeInx), 'Sortstr','descend');
        maxResp(ii,ee) = pks(1);
        maxTmp(ii,ee) = locs(1);
        maxRespSE(ii,ee) = out.broadband.(elecNames{ee}).se(ii,trials.time==locs(1));
    end
end
sumResp = squeeze(sum(elData(:,timeInx,:),2));

% PLOT stimdur versus MAX
figure('Name', [condName ' peak response']);hold on
colors = out.colors;

for ee = 1:length(elecNames)
    subplot(specs.plot.nSubPlots(1),specs.plot.nSubPlots(2),ee);hold on
    for ii = 1:length(whichTrials)
        plot(stimDurations(ii),maxResp(ii,ee),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        e = errorbar(stimDurations(ii),maxResp(ii,ee),maxRespSE(ii,ee), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    %ylim([0 max(maxResp(:)) + 0.1*max(maxResp(:))]);
    ylim([min(maxResp(:)) - 0.1*min(maxResp(:)) max(maxResp(:)) + 0.1*max(maxResp(:))]);

    title(out.titles{ee});
    if ee == 1
        xlabel('stimduration');
        ylabel('peak response level');
        legend(whichTrials, 'Location', 'best')
    end
    set(gca, 'XScale','log', 'YScale','log');
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 1250]);

% PLOT time of max response
figure('Name', [condName ' peak response latency']);hold on

for ee = 1:length(elecNames)
    subplot(specs.plot.nSubPlots(1),specs.plot.nSubPlots(2),ee);hold on
    for ii = 1:length(whichTrials)
        plot(stimDurations(ii),maxTmp(ii,ee),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
    end
	%xlim([0 length(whichTrials)+1]);
    ylim([min(maxTmp(:)) - 0.1*min(maxTmp(:)) max(maxTmp(:)) + 0.1*max(maxTmp(:))]);
    title(out.titles{ee});
    if ee == 1
        xlabel('stim duration');
        ylabel('peak latency (s)');
        legend(whichTrials, 'Location', 'best')
    end

    set(gca, 'XScale','log', 'YScale','log');
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 1250]);

% PLOT SUM
figure('Name', [condName ' sum response']);hold on
colors = out.colors;

for ee = 1:length(elecNames)
    subplot(specs.plot.nSubPlots(1),specs.plot.nSubPlots(2),ee);hold on
    for ii = 1:length(whichTrials)
        plot(stimDurations(ii),sumResp(ii,ee),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        %e = errorbar(CRFrms(ii),maxResp(ii,ee),maxRespSE(ii,ee), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        %e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    %ylim([0 max(maxResp(:)) + 0.1*max(maxResp(:))]);
    ylim([min(sumResp(:)) - 0.1*min(sumResp(:)) max(sumResp(:)) + 0.1*max(sumResp(:))]);

    title(out.titles{ee});
    if ee == 1
        xlabel('stim duration');
        ylabel('sum of response (0-1 sec)');
        legend(whichTrials, 'Location', 'best')
    end
    set(gca, 'XScale','log', 'YScale','log');
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 1250]);

% SAME PLOTs for electrode averages

% time courses
figure('Name', [condName ' timecourses electrode-average']);hold on
for ee = 1:length(eeToAverage)
    elIndex = eeToAverage{ee};
    mnToPlot = mean(elData(:,:,elIndex),3);
    seToPlot = std(elData(:,:,elIndex),0,3);
    subplot(1,3,ee); hold on;
    for ii = 1:length(whichTrials)
        [pks,locs] = findpeaks(mnToPlot(ii,:), trials.time, 'Sortstr','descend');
        pks = pks([1]); 
        locs = locs([1]);
        %c2 = ciplot(mnToPlot(ii,:)-seToPlot(ii,:),mnToPlot(ii,:)+seToPlot(ii,:),trials.time,colors(ii,:), 0.1);
        %c2.Annotation.LegendInformation.IconDisplayStyle = 'off';  
        p1 = plot(trials.time, mnToPlot(ii,:), 'Color', out.colors(ii,:), 'LineWidth',2); 
        p2 = plot(locs,pks, 'Marker', '.', 'MarkerSize', 50, 'Color', colors(ii,:), 'LineStyle', 'none');
        p2.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    end
    l1 = line([0 0], [-10 30],'LineStyle', ':', 'Color', 'k');
    l1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    l2 = line([trials.time(1) trials.time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
    l2.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    xlabel('time(s)');
    ylabel('broadband');    
    set(gca, 'XLim', [-0.2 1.2]);
    tmpYlim = get(gca,'YLim');
    %set(gca, 'YLim', [0 tmpYlim(2)])
    set(gca, 'YLim', [0 20])
    set(gca, 'FontSize', 18);
    title(eeAvTitles{ee});
    set(gcf, 'Position', [150 100 2000 700]);
    if ee == 1
        legend(whichTrials, 'Location', 'best')
    end
end

% max response
figure('Name', [condName ' peak response electrode-average']);hold on

for ee = 1:length(eeToAverage)
    elIndex = eeToAverage{ee};
    maxToPlot = mean(maxResp(:,elIndex),2);
    maxToPlotSE = std(maxResp(:,elIndex),0,2)/sqrt(length(elIndex));
    subplot(1,3,ee); hold on;
    for ii = 1:length(whichTrials)
        plot(stimDurations(ii),maxToPlot(ii),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        e = errorbar(stimDurations(ii),maxToPlot(ii),maxToPlotSE(ii), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    ylim([0 max(maxResp(:)) + 0.1*max(maxResp(:))]);
    title(eeAvTitles{ee});
    if ee == 1
        xlabel('stim duration');
        ylabel('peak response level');
        legend(whichTrials,'Location', 'best')
    end
    set(gca, 'XScale','log', 'YScale','log');    
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 700]);

% time of max response
figure('Name', [condName ' peak response latency electrode-average']);hold on

for ee = 1:length(eeToAverage)
    elIndex = eeToAverage{ee};
    maxToPlot = mean(maxTmp(:,elIndex),2);
    maxToPlotSE = std(maxTmp(:,elIndex),0,2)/sqrt(length(elIndex));
    subplot(1,3,ee); hold on;
    for ii = 1:length(whichTrials)
        plot(stimDurations(ii),maxToPlot(ii),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        e = errorbar(stimDurations(ii),maxToPlot(ii),maxToPlotSE(ii), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    ylim([min(maxTmp(:)) - 0.1*min(maxTmp(:)) max(maxTmp(:)) + 0.1*max(maxTmp(:))]);
    title(eeAvTitles{ee});
    if ee == 1
        xlabel('stim duration');
        ylabel('peak latency (s)');
        legend(whichTrials,'Location', 'best')
    end
    set(gca, 'XScale','log', 'YScale','log')
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 700]);

% sum of response
figure('Name', [condName ' sum response electrode-average']);hold on

for ee = 1:length(eeToAverage)
    elIndex = eeToAverage{ee};
    maxToPlot = mean(sumResp(:,elIndex),2);
    maxToPlotSE = std(sumResp(:,elIndex),0,2)/sqrt(length(elIndex));
    subplot(1,3,ee); hold on;
    for ii = 1:length(whichTrials)
        plot(stimDurations(ii),maxToPlot(ii),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        e = errorbar(stimDurations(ii),maxToPlot(ii),maxToPlotSE(ii), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    ylim([min(sumResp(:)) - 0.1*min(sumResp(:)) max(sumResp(:)) + 0.1*max(sumResp(:))]);
    title(eeAvTitles{ee});
    if ee == 1
        xlabel('stim duration');
        ylabel('sum of response (0-1 sec)');
        legend(whichTrials, 'Location', 'best')
    end
    set(gca, 'XScale','log', 'YScale','log')
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 700]);

% SAVE figures?
if saveFigure
    fg1Nm = [condName '_individualelecs_timecourses'];
    fg2Nm = [condName '_individualelecs_peakresponse'];
    fg3Nm = [condName '_individualelecs_peaklatency'];
    fg4Nm = [condName '_individualelecs_sumresponse'];
    fg5Nm = [condName '_electrodeaverages_timecourses'];
    fg6Nm = [condName '_electrodeaverages_peakresponse'];
    fg7Nm = [condName '_electrodeaverages_peaklatency'];
    fg8Nm = [condName '_electrodeaverages_sumresponse'];
    saveLoc = fullfile(dn_ctrst_RootPath, 'dataFigures');
    saveas(1, [saveLoc filesep fg1Nm], 'epsc');
    saveas(2, [saveLoc filesep fg2Nm], 'epsc');
    saveas(3, [saveLoc filesep fg3Nm], 'epsc');
    saveas(4, [saveLoc filesep fg4Nm], 'epsc');
    saveas(5, [saveLoc filesep fg5Nm], 'epsc');
    saveas(6, [saveLoc filesep fg6Nm], 'epsc');
    saveas(7, [saveLoc filesep fg7Nm], 'epsc');
    saveas(8, [saveLoc filesep fg8Nm], 'epsc');
    %fig = get(1); fig.PaperPosition = fig.PaperPosition/2; print(fig, '-depsc2', [saveLoc filesep fg1Nm '.eps']);
    %printnice(1, 0, saveLoc, fg1Nm)
end
 