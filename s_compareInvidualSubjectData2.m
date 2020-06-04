
conditionsOfInterest = [17];
normalize = 0;
%model = 6;

saveDir = '/Volumes/server/Projects/BAIR/Papers/TemporalDynamicsECoG/figures/temp';
saveFig = 0;

% Plot multiple areas together in one plot
nChans = height(channels);
colors = jet(nChans);

for model = 1:length(results)
    if normalize 
        figName = sprintf('normalized_%s_FullContrast500ms', func2str(results(model).model));
    else
        figName = sprintf('%s_FullContrast500ms', func2str(results(model).model));
    end
    figure('Name', figName);hold on    

    subplot(1,3,1); hold on
    for ii = 1:nChans
        m = squeeze(mean(data2fit(:,conditionsOfInterest,ii),2));
        if normalize, m = m./max(m); end
        if ~isempty(m)
            %plot(t,m,'Color', colors(ii,:), 'LineWidth', 2);        
            if normalize
                ecog_plotSingleTimeCourse(t, m, [], colors(ii,:), 'Data', 'Response normalized', [-0.1 1.1], 'Time');
            else
                ecog_plotSingleTimeCourse(t, m, [], colors(ii,:), 'Data', 'Response', [], 'Time');
            end
        end
    end
    legend(channels.name); 

    subplot(1,3,2); hold on
    for ii = 1:nChans    
        m = squeeze(mean(results(model).pred(:,conditionsOfInterest,ii),2));
        if normalize, m = m/max(m); end
        %plot(t,m,'Color', colors(ii,:), 'LineWidth', 2);  
        if ~isempty(m)
            %plot(t,m,'Color', colors(ii,:), 'LineWidth', 2);        
            if normalize
                ecog_plotSingleTimeCourse(t, m, [], colors(ii,:), [],[], [-0.1 1.1], 'Time');
            else
                ecog_plotSingleTimeCourse(t, m, [], colors(ii,:), [],[], [], 'Time');
            end
        end
    end
    legend(channels.name); xlabel('Time'), title(sprintf('%s Prediction', func2str(results(model).model)));

    if normalize, set(gca, 'Ylim', [-0.1 1.1]);end

    subplot(1,3,3); hold on
    for ii = 1:nChans    
        m = results(model).derived.pred(:,ii);
        if normalize, m = m/max(m); end
        plot(m,'Color', colors(ii,:), 'LineWidth', 2);        
    end
    xlabel('Time'),  title(sprintf('%s derivedPrediction', func2str(results(model).model)))
    legend(channels.name); 

    set(gca, 'Xlim', [0 1200]);
    if normalize, set(gca, 'Ylim', [-0.1 1.1]);end

    set(gcf, 'Position', [600 700 1600 600]);

    set(findall(gcf,'-property','FontSize'),'FontSize',14)

    % Save plot
    if saveFig
        figName = strrep(figName, ' ', '_');
        savePlot(figName, saveDir, 1)
    end
end

% plot saving
function savePlot(figName, saveDir, dataWasAveraged)
    if ~dataWasAveraged
        figDir = fullfile(saveDir, 'individualelectrodes');
    else
        figDir = fullfile(saveDir, 'electrodeaverages');
    end
    if ~exist(figDir, 'dir'), mkdir(figDir), end
    saveas(gcf, fullfile(figDir, figName), 'png'); %close;
    %saveas(gcf, fullfile(figDir, figName), 'fig'); %close;
end