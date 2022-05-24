function tde_plotData(data, channels, t, opts, savePlot, saveStr)

% Generate two plot of selected data:
% 1. each electrode/area in a separate panel, all conditions superimposed; 
% 2. each condition as a separate timecourse, all areas/electrodes superimposed
% Data should be time x conditions x electrodes
%
% 2022 Iris Groen

if ~exist('saveStr', 'var') || isempty(saveStr)
    saveStr = '';
end
if ~exist('savePlot', 'var') || isempty(savePlot)
    savePlot = false;
end

% Plot final selected data
fprintf('[%s] Plotting selected data ... \n',mfilename);

if opts.average_elecs
    figureName1 = sprintf('selecteddata_byarea_electrodeaverages%s', saveStr);
	figureName2 = sprintf('selecteddata_bystimulus_electrodeaverages%s', saveStr);
    FontSz = 20;
    FigSz = [400 200 2000 1200];
else
    figureName1 = sprintf('selecteddata_byarea_individualelecs%s', saveStr);
	figureName2 = sprintf('selecteddata_bystimulus_individualelecs%s', saveStr);
    FontSz = 12;
    FigSz = get(0, 'Screensize');
end

nPlot = size(data,3);
nRow  = ceil(sqrt(nPlot));
nCol  = ceil(sqrt(nPlot));
if nPlot <= (nRow*nCol)-nCol, nRow = nRow-1; end

% Plot each electrode/area in a separate panel, all conditions superimposed
figure('Name', figureName1); 
colors = copper(length(opts.stimnames));
for ii = 1:size(data,3)
    subplot(nRow,nCol,ii); 
    %plot(t,data(:,:,ii), 'LineWidth', 2); colormap(jet);
    
    ecog_plotMultipleTimeCourses(t, data(:,:,ii), [], colors);
    
    if opts.average_elecs
        title(sprintf('%s (n = %d)', ...
            channels.name{ii}, channels.number_of_elecs(ii))); 
    else
        title(sprintf('%s %s %s %s', ...
            channels.benson14_varea{ii}, channels.wang15_mplbl{ii}, channels.subject_name{ii}, channels.name{ii})); 
    end
    if ii == 1
        ylabel('x-fold increase in broadband power');
        xlabel('time (s)');
    end
    yaxlims = get(gca, 'YLim');
    %set(gca, 'YLim', [-0.5 ceil(yaxlims(2)+ (0.1 * yaxlims(2)))]);
    %yaxlims = get(gca, 'YLim');
    %line([0 0], [yaxlims(1) yaxlims(2)], 'Color', 'k', 'LineStyle', ':')
    %line([t(1) t(end)], [0 0],'Color', 'k', 'LineStyle', ':');
    %setsubplotaxes()
    %if ii == 1 && opts.average_elecs
        %legend(opts.stimnames);
        %legend({'ISI 0.017s', 'ISI 0.033s', 'ISI 0.067s', 'ISI 0.133s', 'ISI 0.267s', 'ISI 0.533s'});
    %end
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
    saveas(gcf, fullfile(saveDir, figureName1), 'png'); close;
end

% Plot each condition as a separate timecourse, all areas/electrodes superimposed
figure('Name', figureName2); hold on;
colors = jet(height(channels));
for ii = 1:height(channels), plot(flatten(data(:,:,ii)), 'LineWidth', 2, 'Color', colors(ii,:));end
set(gca, 'xtick', (0:length(opts.stimnames))*size(data,1)+1, 'xgrid', 'on', 'xticklabel', opts.stimnames, 'xticklabelrotation', 45);
axis tight
set(gcf, 'Position', FigSz);
if opts.average_elecs
    legend(channels.name);
end
ylabel('x-fold increase in broadband power');
set(gca, 'FontSize', FontSz);

% save Plot?
if savePlot
    saveDir = fullfile(analysisRootPath, 'figures', 'data');
    if ~exist(saveDir, 'dir'), mkdir(saveDir);end
    fprintf('[%s] Saving figures to %s \n',mfilename, saveDir);
    saveas(gcf, fullfile(saveDir, figureName2), 'png'); close;
end

end