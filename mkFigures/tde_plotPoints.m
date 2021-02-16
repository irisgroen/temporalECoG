function [hp, hc] = tde_plotPoints(m, se, x, type, normalize, line_style)

if ~exist('line_style', 'var') || isempty(line_style)
    line_style = 'none';
end

if normalize
    % Normalize
    se(:,1) = se(:,1)./m(end);
    se(:,2) = se(:,2)./m(end);
    m = m./m(end);
end
switch type
    case 'errbar'
        hp = errorbar(x, m, m-se(:,1), se(:,2)-m, '.', 'Color', 'k', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', line_style, 'CapSize', 0);
        hc = [];
    case 'ci'
        hp = plot(x, m, 'r', 'linewidth', 2);
        hc = ciplot(se(:,1), se(:,2), x, 'r', 0.25);
        set(get(get(hc,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

end