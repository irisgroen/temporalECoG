function tde_plotBars(m, se, x, cmap)

% Plot data as bars with error bars.
%
% 2022 Iris Groen

hb = bar(x,m');

numgroups = size(m,2);
numbars = size(m,1);
groupwidth = min(0.8,numbars/(numbars+1.5));

for ii = 1:numbars
    hb(ii).FaceColor = cmap(ii,:);
end

neg = m-se(:,:,1);
pos = se(:,:,2)-m;
for ii = 1:numbars
    x = (1:numgroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar    
    errorbar(x, m(ii,:), neg(ii,:), pos(ii,:), 'k', 'LineWidth', 2,  'LineStyle', 'none', 'CapSize', 0);
end

end