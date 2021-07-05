% tde_mkFigure 7

modelfuns = {@LINEAR,@LINEAR_RECTF,@LINEAR_RECTF_EXP,@LINEAR_RECTF_EXP_NORM,@LINEAR_RECTF_EXP_NORM_DELAY,...
    @DN, @TTC,@TTCSTIG17, @TTCSTIG19, @HEEGER92, @HEEGER93}; 

xvalmode = 1;
datatype = 'individualelecs';

for ii = 1:length(modelfuns)
    [d(ii)] = tde_loadDataForFigure(modelfuns{ii}, xvalmode, datatype);
end

[results] = tde_evaluateModelFit(d,0);

%% Plot specs

% Subplot positions: % [left bottom width height]
posa = [0.05 0.55 0.85 0.4];
posb = [0.05 0.10 0.85 0.4];

figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

% Panel A: cross-validated R2 for model build up in each area
modelind = [1 2 3 4 5];
R2 = nan(length(modelind),height(d(1).channels));

for ii = 1:size(R2,1)
    %R2(ii,:) = results(modelind(ii)).R2.concat_cond(1,:);
    R2(ii,:) = results(modelind(ii)).R2.concat_all;
end

[~, channels, group_prob] = groupElecsByVisualArea(d(1).channels, 'probabilisticresample');   
[m, se] = averageWithinArea(R2, group_prob, [], 1000);

subplot('position', posa);  cla; hold on
x = 1:height(channels);
cmap = flipud(brewermap(size(R2,1),'RdBu'));
tde_plotBars(m,se,x,cmap);

l = {'DN linear', '+ rectification', '+ exponentiation', '+ normalization', '+ delay'};
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
    %R2(ii,:) = results(modelind(ii)).R2.concat_cond(1,:);
    R2(ii,:) = results(modelind(ii)).R2.concat_all;
end

[~, channels, group_prob] = groupElecsByVisualArea(d(1).channels, 'probabilisticresample');   
[m, se] = averageWithinArea(R2, group_prob, [], 1000);

subplot('position', posb); cla; hold on
x = 1:height(channels);
cmap = brewermap(size(R2,1),'RdYlGn');
tde_plotBars(m,se,x,cmap);

%l = {'LINEAR_RECTF_EXP_NORM_DELAY', 'DN (Zhou 2019)', 'TTC (Horiguchi)', 'TTC (Stigliani 2017)', 'A + S (Stigliani 2019)'};
%l = {'current', 'DN', 'Heeger93', 'TTC 2019', 'TTC 2017', 'TTC'};
l = {'DN (l+r+e+n+d)', 'DN (Zhou19)', 'feedback (Heeger92)', 'feedback (Heeger93)', 'A+S (Stig19)','TTC (Stig17)', 'TTC (Hori09)'};

legend(l, 'location', 'northeast');
legend('boxoff')
set(gca, 'xtick', 1:height(channels), 'xticklabel', channels.name)
set(gca, 'ylim', [0 1], 'xlim', [0 size(m,2)+1]);
box off
ylabel('cross-validated R^{2}', 'Interpreter','tex');
xlabel('visual area');

set(findall(gcf,'-property','FontSize'),'FontSize',20)

% %% temp
% pred = nan(4,2000,height(d(1).channels));
% modelind = [5 6 7 8];
% 
% for ii = 1:size(R2,1)
%     pred(ii,:,:) = results(modelind(ii)).derived.pred;
% end
% 
% [m, se] = averageWithinArea(pred, group_prob, @mean, 1000);
% 
% figure;hold on
% subplot(1,2,1);
% p = plot(m(:,:,1)', 'LineWidth', 2);
% set(p, {'color'}, num2cell(cmap,2));
% %xlim([0 500]);
% legend(l)
%      
% subplot(1,2,2);
% p = plot(m(:,:,1)', 'LineWidth', 2);
% set(p, {'color'}, num2cell(cmap,2));
% xlim([0 400]);
% legend(l)
% 
% 
% figure;hold on;
% for ii = 1:size(m,3)
%     subplot(2,5,ii);
%     p = plot(m(:,:,ii)', 'LineWidth', 2);
%     set(p, {'color'}, num2cell(cmap,2));
% 
%     xlim([0 1000]);
%     if ii == 1, legend(l), end
% end
% 
% figure;hold on
% cmap2 = flipud(brewermap(size(m,3),'YlOrBr'));
% toplot = squeeze(m(1,:,:));
% toplot = toplot./max(toplot);
% p = plot(toplot, 'LineWidth', 2);
% set(p, {'color'}, num2cell(cmap2(1:end,:),2));
% xlim([0 400]);
% legend(channels.name);
% 
% %%
% dat = d(5).data;
% [m, se] = averageWithinArea(dat, group_prob, @mean, 1000);
% 
% figure;hold on
% cmap2 = flipud(brewermap(size(m,3)+2,'YlGnBu'));
% toplot = squeeze(m(:,5,:));
% toplot = toplot./max(toplot);
% p = plot(toplot, 'LineWidth', 2);
% set(p, {'color'}, num2cell(cmap2(1:end-2,:),2));
% xlim([0 400]);
% legend(channels.name);

