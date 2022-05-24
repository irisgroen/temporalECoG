% Figure7-2 : 
% Predicted time courses for deconstructed DN model
%
% 2022 Iris Groen

modelfuns = {@LINEAR,@LINEAR_RECTF,@LINEAR_RECTF_EXP,@LINEAR_RECTF_EXP_NORM,@LINEAR_RECTF_EXP_NORM_DELAY};

xvalmode = 1;
datatype = 'individualelecs';

for ii = 1:length(modelfuns)
    [d(ii)] = tde_loadDataForFigure(modelfuns{ii}, xvalmode, datatype);
end

[results] = tde_evaluateModelFit(d,0);
numboot = 1000;

nModels = length(modelfuns);
cmap = flipud(brewermap(nModels,'Spectral'));

%% OPTION 1
pos = [];

% Subplot positions: % [left bottom width height]
pos(1,:) = [0.05 0.82 0.9 0.15];
pos(2,:) = [0.05 0.62 0.9 0.15];
pos(3,:) = [0.05 0.42 0.9 0.15];
pos(4,:) = [0.05 0.22 0.9 0.15];
pos(5,:) = [0.05 0.02 0.9 0.15];

figure(1); clf; hold on
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
