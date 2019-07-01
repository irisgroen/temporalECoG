%% PLOT SINGLE TRIALS
whichTrials = {'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
whichTrials = {'ONEPULSE-1', 'ONEPULSE-2','ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
whichTrials = {'CRF-1', 'CRF-2','CRF-3', 'CRF-4', 'CRF-5'};

baseline_index = find(~contains(trials.events.trial_name, 'PRF'));
                
for ee = 1:length(elecNames)
    elecInd = ecog_matchChannels(elecNames{ee}, trials);
    baseline = mean(mean(squeeze(trials.broadband(elecInd,trials.time<0,baseline_index)),1),2);
    figure('Name', out.titles{ee});hold on;
    for ii = 1:length(whichTrials)
        subplot(6,1,ii);
        stimInd = contains(trials.events.trial_name, whichTrials{ii});
        trialsToPlot = squeeze(trials.broadband(elecInd,:,stimInd));
        trialsToPlot = (trialsToPlot - baseline) ./ baseline;
        plot(trials.time, trialsToPlot, 'LineWidth', 2);
        title(whichTrials{ii});
        set(gca, 'XLim', [-0.2 1.2]);
        tmpYlim = get(gca,'YLim');
        %set(gca, 'YLim', [0 tmpYlim(2)])
        set(gca, 'YLim', [0 50])
        set(gca, 'FontSize', 18);
        xlabel('time(s)');
        ylabel('broadband');
    end
    set(gcf, 'Position', [200 0 750 1300]);
end