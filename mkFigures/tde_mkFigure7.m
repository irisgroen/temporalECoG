% tde_mkFigure 7

modelfuns = {@LINEAR,@LINEAR_RECTF,@LINEAR_RECTF_EXP,@LINEAR_RECTF_EXP_NORM,@LINEAR_RECTF_EXP_NORM_DELAY,...
    @DN, @TTC,@TTCSTIG17, @TTCSTIG19, @HEEGER92, @HEEGER93}; 

xvalmode = 1;
datatype = 'individualelecs';

for ii = 1:length(modelfuns)
    [d(ii)] = tde_loadDataForFigure(modelfuns{ii}, xvalmode, datatype);
end

[results] = tde_evaluateModelFit(d,0);
numboot = 10000;

%% Plot specs

% Subplot positions: % [left bottom width height]
posa = [0.05 0.55 0.95 0.4];
posb = [0.05 0.10 0.95 0.4];

figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

% Panel A: cross-validated R2 for model build up in each area
modelind = [1 2 3 4 5];
R2 = nan(length(modelind),height(d(1).channels));

for ii = 1:size(R2,1)
    R2(ii,:) = results(modelind(ii)).R2.concat_all;
end

[~, channels, group_prob] = groupElecsByVisualArea(d(1).channels, 'probabilisticresample');   
[m, se] = averageWithinArea(R2, group_prob, [], numboot);

subplot('position', posa);  cla; hold on
x = 1:height(channels);
cmap = flipud(brewermap(size(R2,1),'RdBu'));
tde_plotBars(m,se,x,cmap);

l = {'linear', '+ rectification', '+ exponentiation', '+ normalization', '+ delay'};
legend(l, 'location', 'northeast');
legend('boxoff')
set(gca, 'xtick', 1:height(channels), 'xticklabel', channels.name)
set(gca, 'ylim', [0 1], 'xlim', [0 size(m,2)+1]);
box off
ylabel('cross-validated R^{2}', 'Interpreter','tex'); 

set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% Panel B: cross-validated R2 for other models vs DN in all areas

modelind = [5 6 10 11 9 8 7];
R2 = nan(length(modelind),height(d(1).channels));

for ii = 1:size(R2,1)
    R2(ii,:) = results(modelind(ii)).R2.concat_all;
end

[~, channels, group_prob] = groupElecsByVisualArea(d(1).channels, 'probabilisticresample');   
[m, se] = averageWithinArea(R2, group_prob, [], numboot);

subplot('position', posb); cla; hold on
x = 1:height(channels);
cmap = brewermap(size(R2,1),'RdYlGn');
tde_plotBars(m,se,x,cmap);

l = {'DN (flexible IRF)', 'DN (constrained IRF)', 'feedback (Heeger92)', 'feedback (Heeger93)', 'A+S (Stig19)','TTC (Stig17)', 'TTC (Hori09)'};

legend(l, 'location', 'northeast');
legend('boxoff')
set(gca, 'xtick', 1:height(channels), 'xticklabel', channels.name)
set(gca, 'ylim', [0 1], 'xlim', [0 size(m,2)+1]);
box off
ylabel('cross-validated R^{2}', 'Interpreter','tex');
xlabel('visual area');

set(findall(gcf,'-property','FontSize'),'FontSize',24)


