function tde_plotData(data, channels, t, opts, savePlot)

% Data should be time x trials x channels

if ~exist('savePlot', 'var') || isempty(savePlot)
    savePlot = 0;
end

% Plot final selected data
fprintf('[%s] Plotting selected data ... \n',mfilename);

if opts.average_elecs
    figureName = sprintf('selecteddata_electrodeaverages_%s', opts.stimnames{1}(1:end-2));
    FontSz = 20;
    FigSz = [400 200 2000 1200];
else
    figureName = sprintf('selecteddata_individualelecs_%s', opts.stimnames{1}(1:end-2));
    FontSz = 12;
    FigSz = get(0, 'Screensize');
end

nPlot = size(data,3);
nRow  = ceil(sqrt(nPlot));
nCol  = ceil(sqrt(nPlot));
if nPlot <= (nRow*nCol)-nCol, nRow = nRow-1; end

figure('Name', figureName); 
colors = copper(length(opts.stimnames));
for ii = 1:size(data,3)
    subplot(nRow,nCol,ii); 
    %plot(t,data(:,:,ii), 'LineWidth', 2); colormap(jet);
    for ss = 1:size(data,2)
        ecog_plotSingleTimeCourse(t, data(:,ss,ii), [], colors(ss,:));
    end

    if opts.average_elecs
        title(sprintf('%s (n = %d)', ...
            channels.name{ii}, channels.number_of_elecs{ii})); 
    else
        title(sprintf('%s %s %s %s', ...
            channels.bensonarea{ii}, channels.wangarea{ii}, channels.subject_name{ii}, channels.name{ii})); 
    end
    if ii == 1
        ylabel('x-fold increase in broadband power');
        xlabel('time (s)');
    end
    yaxlims = get(gca, 'YLim');
    set(gca, 'YLim', [-0.5 ceil(yaxlims(2)+ (0.1 * yaxlims(2)))]);
    %yaxlims = get(gca, 'YLim');
    %line([0 0], [yaxlims(1) yaxlims(2)], 'Color', 'k', 'LineStyle', ':')
    %line([t(1) t(end)], [0 0],'Color', 'k', 'LineStyle', ':');
    %setsubplotaxes()
    if ii == 1 && opts.average_elecs
        %legend(opts.stimnames);
        legend({'ISI 0.017s', 'ISI 0.033s', 'ISI 0.067s', 'ISI 0.133s', 'ISI 0.267s', 'ISI 0.533s'});
    end
    set(gca, 'FontSize', FontSz);
    
    if ~opts.average_elecs
        set(gca, 'XTickLabel', [], 'YTickLabel', []);
    end
end
set(gcf, 'Position', FigSz);

% save Plot?
if savePlot
    saveDir = fullfile(analysisRootPath, 'figures', 'data');
    if ~exist(saveDir, 'dir'), mkdir(saveDir);end
    fprintf('[%s] Saving figures to %s \n',mfilename, saveDir);
    saveas(gcf, fullfile(saveDir, figureName), 'png');
end
end