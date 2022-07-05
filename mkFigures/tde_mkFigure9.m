% Figure 9 : 
% Delayed normalization is especially important for explaining contrast-related dynamics
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

%% Panel A (here Figure 1)

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


%% Panel B (here figure 2)

numboot = 1000;

pos = [];

% Subplot positions: % [left bottom width height]
pos(1,:) = [0.05 0.82 0.9 0.15];
pos(2,:) = [0.05 0.62 0.9 0.15];
pos(3,:) = [0.05 0.42 0.9 0.15];
pos(4,:) = [0.05 0.22 0.9 0.15];
pos(5,:) = [0.05 0.02 0.9 0.15];

figure(2); clf; hold on
set(gcf, 'position',  get(0, 'screensize'));

[~, ~, group_prob] = groupElecsByVisualArea(d(1).channels, 'probabilisticresample', {'V1'});   

conditionsOfInterest = {'ONEPULSE', 'TWOPULSE', 'CRF'}; 
timepointsOfInterest = [-0.10 1];

% Look up corresponding indices in the data
stim_idx = [];
for ii = 1:length(conditionsOfInterest), stim_idx = [stim_idx; find(contains(d(1).stim_info.name, conditionsOfInterest{ii}))]; end
t_idx    = d(1).t>timepointsOfInterest(1) & d(1).t<=timepointsOfInterest(2);

for ii = 1:5

    subplot('position', pos(ii,:)); cla; hold on
    D = d(ii);

    stim = D.stim(t_idx,stim_idx);
    [data, data_se] = averageWithinArea(D.data(t_idx,stim_idx,:), group_prob, @mean, numboot);
    [pred, pred_se] = averageWithinArea(D.pred(t_idx,stim_idx,:), group_prob, @mean, numboot);

    s = flatten(stim);
    plot(s,'color', [0.5 0.5 0.5], 'lineWidth', 1); 
    
    dat = flatten(data);
    plot(smooth(dat/max(dat,[],2),15),'k', 'LineWidth',2); 
    
    prd = flatten(pred);
    plot(smooth(prd/max(prd,[],2),15),'Color',cmap(ii,:), 'LineWidth',2); 

    set(gca, 'XTick',1:length(find(t_idx)):length(stim_idx)*length(find(t_idx)), 'XTickLabel', []);
    axis tight
    set(gca, 'YLim', [-0.3 1]);
    set(gca, 'ytick', [0 1]);
end

set(findall(gcf,'-property','FontSize'),'FontSize',24)
print('-painters','-dsvg','/Users/iiagroen/surfdrive/BAIR/Papers/TemporalDynamicsECoG/mkFigures/revision/ExtendedData_Figure7_2.svg')

