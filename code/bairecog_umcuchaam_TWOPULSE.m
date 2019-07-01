saveFigure = 1;

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

eeToAverage = {[1:6],[6:8],[1:8]}; 
eeAvTitles = {'all V1 electrodes', 'all V2/V3 electrodes', 'all V1/V2/V3 electrodes'};
elecNames = trials.channels.name;

%% TWOPULSE %%
close all;

% PLOT time course for each condition
specs = [];
specs.dataTypes          = {'broadband'};
specs.smoothingLevelInMs = [];

specs.plot.colorMap      = 'parula';
specs.plot.nSubPlots     = [2 4];
specs.plot.addEccToTitle = 'yes';
specs.plot.showMax       = 'yes';

whichTrials = {'ONEPULSE-5', 'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
[out] = ecog_plotTimecourses(trials, elecNames, whichTrials,specs);

condName = whichTrials{2}(1:8);
stimDurations = unique(stimData.stimulus.ISI);

%%
close all;

colors = out.colors;
timeInd = trials.time>0 & trials.time<1;

% Detect peaks and location of peaks; make plots for inspection
allPks = []; allLocs = []; allMn = [];
for ee = 1:length(elecNames)
    mnToPlot = out.broadband.(elecNames{ee}).mn;
    seToPlot = out.broadband.(elecNames{ee}).se;
    figure('Name', out.titles{ee});hold on;

    for ii = 1:length(whichTrials)
        [pks,locs] = findpeaks(mnToPlot(ii,timeInd), trials.time(timeInd), ...
            'MinPeakDistance',stimDurations(ii)+mode(stimData.stimulus.duration)/2, 'Sortstr','descend');
        pks = pks([1 2]); % take top 2 peaks
        locs = locs([1 2]);
        if locs(2) < locs(1)
            pks = pks([2 1]);
            locs = locs([2 1]);
        end       
        if ii > 1
            subplot(3,3,ii+2);
        else
            subplot(3,3,ii);
        end
        hold on;
        p1 = plot(trials.time, mnToPlot(ii,:), 'Color', colors(ii,:), 'LineWidth',2);
        c2 = ciplot(mnToPlot(ii,:)-seToPlot(ii,:),mnToPlot(ii,:)+seToPlot(ii,:),trials.time,colors(ii,:), 0.25);
        c2.Annotation.LegendInformation.IconDisplayStyle = 'off';  
        p2 = plot(locs,pks, 'Marker', '.', 'MarkerSize', 50, 'Color', colors(ii,:), 'LineStyle', 'none');
        p2.Annotation.LegendInformation.IconDisplayStyle = 'off';  
        title(whichTrials{ii});
        l1 = line([0 0], [-10 30],'LineStyle', ':', 'Color', 'k');
        l1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
        l2 = line([trials.time(1) trials.time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
        l2.Annotation.LegendInformation.IconDisplayStyle = 'off';  
        xlabel('time(s)');
        ylabel('broadband');
        
        set(gca, 'XLim', [-0.2 1.2]);
        tmpYlim = get(gca,'YLim');
        %set(gca, 'YLim', [0 tmpYlim(2)])
        set(gca, 'YLim', [-2 30])
        set(gca, 'FontSize', 18);
        
        % save pks and locs
        allPks(ee,ii,:) = pks;
        allLocs(ee,ii,:) = locs;
        
    end
    allMn(ee,:,:) = mnToPlot;
    set(gcf, 'Position', [150 100 2000 1250]);
end

% Make same plot for AVERAGE across elecs
figure('Name', [condName ' time courses electrode-average']);hold on

for ee = 1:length(eeToAverage)
    elIndex = eeToAverage{ee};
 
    mnToPlot = squeeze(mean(allMn(elIndex,:,:),1));
    seToPlot = squeeze(std(allMn(elIndex,:,:),0,1))/sqrt(length(elIndex));
    
    subplot(1,3,ee);hold on;
    for ii = 1:length(whichTrials)
        [pks,locs] = findpeaks(mnToPlot(ii,timeInd), trials.time(timeInd), ...
            'MinPeakDistance',stimDurations(ii)+mode(stimData.stimulus.duration)/2, 'Sortstr','descend');
        pks = pks([1 2]); % take top 2 peaks
        locs = locs([1 2]);
        c2 = ciplot(mnToPlot(ii,:)-seToPlot(ii,:),mnToPlot(ii,:)+seToPlot(ii,:),trials.time,colors(ii,:), 0.25);
        c2.Annotation.LegendInformation.IconDisplayStyle = 'off';  
        p1 = plot(trials.time, mnToPlot(ii,:), 'Color', colors(ii,:), 'LineWidth',2);
        p2 = plot(locs,pks, 'Marker', '.', 'MarkerSize', 50, 'Color', colors(ii,:), 'LineStyle', 'none');
        p2.Annotation.LegendInformation.IconDisplayStyle = 'off';  
        l1 = line([0 0], [-10 30],'LineStyle', ':', 'Color', 'k');
        l1.Annotation.LegendInformation.IconDisplayStyle = 'off';  
        l2 = line([trials.time(1) trials.time(end)], [0 0],'LineStyle', ':', 'Color', 'k');
        l2.Annotation.LegendInformation.IconDisplayStyle = 'off';  
        xlabel('time(s)');
        ylabel('broadband');

        set(gca, 'XLim', [-0.2 1.2]);
        tmpYlim = get(gca,'YLim');
        %set(gca, 'YLim', [0 tmpYlim(2)])
        set(gca, 'YLim', [-2 30])
        set(gca, 'FontSize', 18);
    end   
    title(eeAvTitles{ee});

    set(gcf, 'Position', [150 100 2000 700]);
    if ee == 1
        legend(whichTrials, 'Location', 'NorthWest')
    end
end

% SAVE figures?

if saveFigure
	saveLoc = fullfile(dn_ctrst_RootPath, 'dataFigures');
    for ee = 1:length(elecNames)
        fgNm = [condName '_elec' elecNames{ee} '_timecourses'];
        saveas(ee, [saveLoc filesep fgNm], 'epsc');
    end
    fg1Nm = [condName '_electrodeaverages_timecourses'];
    fg2Nm = [condName '_electrodeaverages_interpeakdistance'];   
    saveas(ee+1, [saveLoc filesep fg1Nm], 'epsc');
    %saveas(ee+2, [saveLoc filesep fg2Nm], 'epsc');
    %fig = get(1); fig.PaperPosition = fig.PaperPosition/2; print(fig, '-depsc2', [saveLoc filesep fg1Nm '.eps']);
    %printnice(1, 0, saveLoc, fg1Nm)
end

%%  more plots
% PLOT first peak latency + max
% PLOT second peak latency + max
% PLOT interpeak distance as function of ISI
legendToPlot = {'first peak', 'second peak'};

close all;

% PLOT stimdur versus peak response level
figure('Name', [condName ' peak response']); clf; hold on
colors = out.colors;

for ee = 1:length(elecNames)
    subplot(specs.plot.nSubPlots(1),specs.plot.nSubPlots(2),ee);hold on
    for ii = 1:length(whichTrials)
        p1 = plot(stimDurations(ii),allPks(ee,ii,1),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        p2 = plot(stimDurations(ii),allPks(ee,ii,2),'Color', colors(ii,:), 'Marker', 'o', 'MarkerSize', 15,'LineStyle', 'none');
        %p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        %e = errorbar(stimDurations(ii),maxResp(ii,ee),maxRespSE(ii,ee), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        %e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    %ylim([0 max(maxResp(:)) + 0.1*max(maxResp(:))]);
    ylim([min(allPks(ee,:)) - 0.1*min(allPks(ee,:)) max(allPks(ee,:)) + 0.1*max(allPks(ee,:))]);

    title(out.titles{ee});
    if ee == 1
        xlabel('stimulus ISI');
        ylabel('peak response level');
        legend(legendToPlot, 'Location', 'best')
    end
    %set(gca, 'XScale','log', 'YScale','log')
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 1250]);

% PLOT stimulus ISI versus peak latency
figure('Name', [condName ' peak response latency']); clf; hold on
colors = out.colors;

for ee = 1:length(elecNames)
    subplot(specs.plot.nSubPlots(1),specs.plot.nSubPlots(2),ee);hold on
    for ii = 1:length(whichTrials)
        p1 = plot(stimDurations(ii),allLocs(ee,ii,1),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        p2 = plot(stimDurations(ii),allLocs(ee,ii,2),'Color', colors(ii,:), 'Marker', 'o', 'MarkerSize', 15,'LineStyle', 'none');
        %p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        %e = errorbar(stimDurations(ii),maxResp(ii,ee),maxRespSE(ii,ee), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        %e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    %ylim([0 max(maxResp(:)) + 0.1*max(maxResp(:))]);
    ylim([min(allLocs(ee,:)) - 0.1*min(allLocs(ee,:)) max(allLocs(ee,:)) + 0.1*max(allLocs(ee,:))]);

    title(out.titles{ee});
    if ee == 1
        xlabel('stimulus ISI');
        ylabel('peak latency (s)');
        legend(legendToPlot, 'Location', 'best')
    end
    %set(gca, 'XScale','log', 'YScale','log')
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 1250]);

% PLOT stimulus ISI versus SUM
timeInx = trials.time>0 & trials.time<1; 
elData = [];
for ee = 1:length(elecNames)
    elData(:,:,ee) = out.broadband.(elecNames{ee}).mn;
end
sumResp = squeeze(sum(elData(:,timeInx,:),2));
figure('Name', [condName ' sum response']);hold on

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
        xlabel('stimulus ISI');
        ylabel('sum of response (0-1 sec)');
        legend(whichTrials, 'Location', 'best')
    end
    %set(gca, 'XScale','log', 'YScale','log')
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 1250]);

% PLOT stimulus ISI versus inter-peak difference in LEVEL
figure('Name', [condName ' inter-peak latency']); clf; hold on
colors = out.colors;
PksToPlot = allPks(:,:,2)-allPks(:,:,1);

for ee = 1:length(elecNames)
    subplot(specs.plot.nSubPlots(1),specs.plot.nSubPlots(2),ee);hold on
    for ii = 1:length(whichTrials)
        p1 = plot(stimDurations(ii),PksToPlot(ee,ii),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        %p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        %e = errorbar(stimDurations(ii),maxResp(ii,ee),maxRespSE(ii,ee), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        %e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    %ylim([0 max(maxResp(:)) + 0.1*max(maxResp(:))]);
    ylim([min(PksToPlot(ee,:)) - 0.1*min(PksToPlot(ee,:)) max(PksToPlot(ee,:)) + 0.1*max(PksToPlot(ee,:))]);

    title(out.titles{ee});
    if ee == 1
        xlabel('stimulus ISI');
        ylabel('inter-peak response difference');
        legend(whichTrials, 'Location', 'best')
    end
    %set(gca, 'XScale','log', 'YScale','log')
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 1250]);

% PLOT stimulus ISI versus inter-peak difference in LATENCY
figure('Name', [condName ' inter-peak latency']); clf; hold on
colors = out.colors;
LocsToPlot = allLocs(:,:,2)-allLocs(:,:,1);

for ee = 1:length(elecNames)
    subplot(specs.plot.nSubPlots(1),specs.plot.nSubPlots(2),ee);hold on
    for ii = 1:length(whichTrials)
        p1 = plot(stimDurations(ii),LocsToPlot(ee,ii),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        %p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
        %e = errorbar(stimDurations(ii),maxResp(ii,ee),maxRespSE(ii,ee), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        %e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    %ylim([0 max(maxResp(:)) + 0.1*max(maxResp(:))]);
    ylim([min(LocsToPlot(ee,:)) - 0.1*min(LocsToPlot(ee,:)) max(LocsToPlot(ee,:)) + 0.1*max(LocsToPlot(ee,:))]);

    title(out.titles{ee});
    if ee == 1
        xlabel('stimulus ISI');
        ylabel('inter-peak latency (s)');
        legend(whichTrials, 'Location', 'best')
    end
    %set(gca, 'XScale','log', 'YScale','log')
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 1250]);


% SAME PLOTs for electrode averages

% peak response
figure('Name', [condName ' peak response electrode-average']);hold on

for ee = 1:length(eeToAverage)
    elIndex = eeToAverage{ee};
         
    maxToPlot = squeeze(mean(allPks(elIndex,:,:),1));
    maxToPlotSE = squeeze(std(allPks(elIndex,:,:),0,1))/sqrt(length(elIndex));
    subplot(1,3,ee); hold on;
    for ii = 1:length(whichTrials)
        p1 = plot(stimDurations(ii),maxToPlot(ii,1),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        p2 = plot(stimDurations(ii),maxToPlot(ii,2),'Color', colors(ii,:), 'Marker', 'o', 'MarkerSize', 15,'LineStyle', 'none');       
        e1 = errorbar(stimDurations(ii),maxToPlot(ii,1),maxToPlotSE(ii,1), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        e2 = errorbar(stimDurations(ii),maxToPlot(ii,2),maxToPlotSE(ii,2), 'LineWidth', 1,'LineStyle', 'none', 'Color', colors(ii,:));
        e1.Annotation.LegendInformation.IconDisplayStyle = 'off';
        e2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    ylim([0 max(allPks(:)) + 0.1*max(allPks(:))]);
    title(eeAvTitles{ee});
    if ee == 1
        xlabel('stim ISI');
        ylabel('response level');
        legend(legendToPlot,'Location', 'best')
    end
    %set(gca, 'XScale','log', 'YScale','log')
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 700]);


% peak latency
figure('Name', [condName ' peak response electrode-average']);hold on

for ee = 1:length(eeToAverage)
    elIndex = eeToAverage{ee};
         
    maxToPlot = squeeze(mean(allLocs(elIndex,:,:),1));
    maxToPlotSE = squeeze(std(allLocs(elIndex,:,:),0,1))/sqrt(length(elIndex));
    subplot(1,3,ee); hold on;
    for ii = 1:length(whichTrials)
        p1 = plot(stimDurations(ii),maxToPlot(ii,1),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        p2 = plot(stimDurations(ii),maxToPlot(ii,2),'Color', colors(ii,:), 'Marker', 'o', 'MarkerSize', 15,'LineStyle', 'none');       
        e1 = errorbar(stimDurations(ii),maxToPlot(ii,1),maxToPlotSE(ii,1), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        e2 = errorbar(stimDurations(ii),maxToPlot(ii,2),maxToPlotSE(ii,2), 'LineWidth', 1,'LineStyle', 'none', 'Color', colors(ii,:));
        e1.Annotation.LegendInformation.IconDisplayStyle = 'off';
        e2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    ylim([0 max(allLocs(:)) + 0.1*max(allLocs(:))]);
    title(eeAvTitles{ee});
    if ee == 1
        xlabel('stim ISI');
        ylabel('response latency');
        legend(legendToPlot,'Location', 'best')
    end
    %set(gca, 'XScale','log', 'YScale','log')
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
        xlabel('stim ISI');
        ylabel('sum of response (0-1 sec)');
        legend(whichTrials, 'Location', 'best')
    end
    %set(gca, 'XScale','log', 'YScale','log')
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 700]);

% PLOT stimulus ISI versus inter-peak difference in LEVEL
figure('Name', [condName ' inter-peak response level electrode-average']); hold on
colors = out.colors;
PksToPlot = allPks(:,:,2)-allPks(:,:,1);

for ee = 1:length(eeToAverage)
    elIndex = eeToAverage{ee};
    maxToPlot = mean(PksToPlot(elIndex,:),1);
    maxToPlotSE = std(PksToPlot(elIndex,:),0,1)/sqrt(length(elIndex));
    subplot(1,3,ee); hold on;
    for ii = 1:length(whichTrials)
        plot(stimDurations(ii),maxToPlot(ii),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        e = errorbar(stimDurations(ii),maxToPlot(ii),maxToPlotSE(ii), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    ylim([min(PksToPlot(:)) - 0.1*min(PksToPlot(:)) max(PksToPlot(:)) + 0.1*max(PksToPlot(:))]);
    title(eeAvTitles{ee});
    if ee == 1
        xlabel('stim ISI');
        ylabel('inter-peak difference in response level');
        legend(whichTrials, 'Location', 'best')
    end
    %set(gca, 'XScale','log', 'YScale','log')
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 700]);

% PLOT stimulus ISI versus inter-peak difference in LATENCY
figure('Name', [condName ' inter-peak latency electrode-average']); hold on
colors = out.colors;
PksToPlot = allLocs(:,:,2)-allLocs(:,:,1);

for ee = 1:length(eeToAverage)
    elIndex = eeToAverage{ee};
    maxToPlot = mean(PksToPlot(elIndex,:),1);
    maxToPlotSE = std(PksToPlot(elIndex,:),0,1)/sqrt(length(elIndex));
    subplot(1,3,ee); hold on;
    for ii = 1:length(whichTrials)
        plot(stimDurations(ii),maxToPlot(ii),'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 50,'LineStyle', 'none');
        e = errorbar(stimDurations(ii),maxToPlot(ii),maxToPlotSE(ii), 'LineWidth', 2,'LineStyle', 'none', 'Color', colors(ii,:));
        e.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
	%xlim([0 length(whichTrials)+1]);
    ylim([min(PksToPlot(:)) - 0.1*min(PksToPlot(:)) max(PksToPlot(:)) + 0.1*max(PksToPlot(:))]);
    title(eeAvTitles{ee});
    if ee == 1
        xlabel('stim ISI');
        ylabel('inter-peak difference in response latency');
        legend(whichTrials, 'Location', 'best')
    end
    %set(gca, 'XScale','log', 'YScale','log')
    set(gca, 'FontSize', 18);
end
set(gcf, 'Position', [150 100 2000 700]);

if saveFigure
	saveLoc = fullfile(dn_ctrst_RootPath, 'dataFigures');
    fg1Nm = [condName '_individualelecs_peakresponse'];
    fg2Nm = [condName '_individualelecs_peaklatency']; 
    fg3Nm = [condName '_individualelecs_sumresponse'];
    fg4Nm = [condName '_individualelecs_interpeakresponsediff'];
    fg5Nm = [condName '_individualelecs_interpeaklatencydiff'];
    fg6Nm = [condName '_electrodeaverages_peakresponse'];
    fg7Nm = [condName '_electrodeaverages_peaklatency']; 
    fg8Nm = [condName '_electrodeaverages_sumresponse'];
    fg9Nm = [condName '_electrodeaverages_interpeakresponsediff'];
    fg10Nm = [condName '_electrodeaverages_interpeaklatencydiff'];
    saveas(1, [saveLoc filesep fg1Nm], 'epsc');
    saveas(2, [saveLoc filesep fg2Nm], 'epsc');
    saveas(3, [saveLoc filesep fg3Nm], 'epsc');
    saveas(4, [saveLoc filesep fg4Nm], 'epsc');
    saveas(5, [saveLoc filesep fg5Nm], 'epsc');
    saveas(6, [saveLoc filesep fg6Nm], 'epsc');
    saveas(7, [saveLoc filesep fg7Nm], 'epsc');
    saveas(8, [saveLoc filesep fg8Nm], 'epsc');
    saveas(9, [saveLoc filesep fg9Nm], 'epsc');
    saveas(10, [saveLoc filesep fg10Nm], 'epsc');
end

