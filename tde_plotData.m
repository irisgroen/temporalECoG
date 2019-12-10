function tde_plotData(data, channels, t, opts)

% Data should be time x trials x channels

% Plot final selected data
figure('Name', 'Selected data'); 
sz = ceil(sqrt(size(data,3)));
for ii = 1:size(data,3)
    subplot(sz,sz,ii); plot(t,data(:,:,ii), 'LineWidth', 2);
    if opts.average_elecs
        title(sprintf('%s (n = %d)', ...
            channels.name{ii}, channels.number_of_elecs{ii})); 
    else
        title(sprintf('%s %s %s %s', ...
            channels.bensonarea{ii}, channels.wangarea{ii}, channels.subject_name{ii}, channels.name{ii})); 
    end
    yaxlims = get(gca, 'YLim');
    set(gca, 'YLim', [-0.5 ceil(yaxlims(2)+ (0.1 * yaxlims(2)))]);
    yaxlims = get(gca, 'YLim');
    line([0 0], [yaxlims(1) yaxlims(2)], 'Color', 'k', 'LineStyle', ':')
    line([t(1) t(end)], [0 0],'Color', 'k', 'LineStyle', ':');
end
set(gcf, 'Position', [400 200 2000 1200]);
end