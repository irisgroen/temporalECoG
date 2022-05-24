% Figure7-1 : 
% Variance explained by DN model, separated by stimulus condition
%
% 2022 Iris Groen

modelfuns = {@LINEAR,@LINEAR_RECTF,@LINEAR_RECTF_EXP,@LINEAR_RECTF_EXP_NORM,@LINEAR_RECTF_EXP_NORM_DELAY,@DN};

xvalmode = 1;
datatype = 'individualelecs';

for ii = 1:length(modelfuns)
    [d(ii)] = tde_loadDataForFigure(modelfuns{ii}, xvalmode, datatype);
end

[results] = tde_evaluateModelFit(d,0);
numboot = 10000;

nModels = length(modelfuns);
cmap = flipud(brewermap(nModels,'Spectral'));

%% Figure 

% Subplot positions: % [left bottom width height]
pos(1,:) = [0.05 0.7 0.95 0.25];
pos(2,:) = [0.05 0.38 0.95 0.25];
pos(3,:) = [0.05 0.07 0.95 0.25];

figure(1); clf; hold on
set(gcf, 'position',  get(0, 'screensize'));

% Cross-validated R2 for model build up in each area
modelind = [1 2 3 4 5 6];
R2 = nan(length(modelind),height(d(1).channels));
conditionNames = {'Contrast', 'Duration', 'Repetition'};
conditionIndex = [2 3 1];    

for jj = 1:3
    
    for ii = 1:size(R2,1)
        R2(ii,:) = results(modelind(ii)).R2.concat_cond(conditionIndex(jj),:);
    end

    [~, channels, group_prob] = groupElecsByVisualArea(d(1).channels, 'probabilisticresample');   
    [m, se] = averageWithinArea(R2, group_prob, [], numboot);

    subplot('position', pos(jj,:));  cla; hold on
    x = 1:height(channels);
    tde_plotBars(m,se,x,cmap);
    
    if jj == 1
        l = {'Linear', '+ Rectification', '+ Exponentiation', '+ Normalization', '+ Delay', '+ Delay (fixed IRF)'};
        legend(l, 'location', 'northeast');
        legend('boxoff')
    end
    set(gca, 'xtick', 1:height(channels));
    if jj == 3, set(gca, 'xticklabel', channels.name), else,  set(gca, 'xticklabel', []), end
    set(gca, 'ylim', [0 1], 'xlim', [0 size(m,2)+1]);
    box off
    ylabel('cross-validated R^{2}', 'Interpreter','tex'); 
    title(conditionNames(conditionIndex(jj)));
end

set(findall(gcf,'-property','FontSize'),'FontSize',24)



