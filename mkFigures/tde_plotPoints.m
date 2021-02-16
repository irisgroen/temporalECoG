function tde_plotPoints(m, se, x, type, normalize)

if normalize
    % Normalize
    se(:,1) = se(:,1)./m(end);
    se(:,2) = se(:,2)./m(end);
    m = m./m(end);
end
switch type
    case 'data'
        errorbar(x, m, m-se(:,1), se(:,2)-m, '.', 'Color', 'k', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0);
    case 'model'
        plot(x, m, 'r', 'linewidth', 2);
        ch = ciplot(se(:,1), se(:,2), x, 'r', 0.25);
        set(get(get(ch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
end

end