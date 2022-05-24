function [hp, hc] = tde_plotPoints(m, se, x, type, normalize, line_style, msize, color)

% Plot data as points with error bars (type = 'errbar'), or as a line with
% shaded confidence regions (type = 'ci').
%
% 2022 Iris Groen

if ~exist('color', 'var') || isempty(color)
    color = 'k';
end

if ~exist('line_style', 'var') || isempty(line_style)
    line_style = 'none';
end

if ~exist('msize', 'var') || isempty(msize)
    msize = 30;
end

if normalize
    % Normalize
    if ~isempty(se)
        se(:,1) = se(:,1)./m(end);
        se(:,2) = se(:,2)./m(end);
    end
    m = m./m(end);
end
switch type
    case 'errbar'
        hp = errorbar(x, m, m-se(:,1), se(:,2)-m, '.', 'Color', color, 'MarkerSize', msize, 'LineWidth', 2, 'LineStyle', line_style, 'CapSize', 0);
        hc = [];
    case 'ci'
        hp = plot(x, m, 'color', color, 'linewidth', 2);
        if ~isempty(se)
            hc = ciplot(se(:,1), se(:,2), x, color, 0.25);
            set(get(get(hc,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
end

end