% Figure 12 :
% Temporal summation windows vary across individual participants
%
% 2022 Iris Groen

% Load data and fits
modelfun = @DN;
xvalmode = 1;
datatype = 'individualelecs';
D = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% Compute derived parameters ysing the model
[results] = tde_evaluateModelFit(D);

% Compute derived parameters using the data
t = D.t;
stim_info = D.stim_info;
COI = {'ONEPULSE'};
stim_idx  = contains(stim_info.name, COI);

[derivedPrm, derivedPrmNames] = tde_computeDerivedParamsData(D.data(:,stim_idx,:),t); 

% Get list of areas 
[~, channels, group_prob] = groupElecsByVisualArea(D.channels, 'probabilisticresample');   

%% Extended data: individual subjects
figure(2); clf
ScrSz = get(0, 'screensize');
set(gcf, 'position',  [ScrSz(1) ScrSz(2) ScrSz(3) ScrSz(4)/2]);
hold on;
x = 1:height(channels);
fun = @median;
nboot = 1000;

ylabels = derivedPrmNames;

xlabel('Visual area');

subjectNames = unique(D.channels.subject_name);

% Recode subjectnames (use for plot titles only)
subjectNames_recoded = subjectNames;
idx = contains(subjectNames, 'p02'); subjectNames_recoded{idx} = 'Patient 2';
idx = contains(subjectNames, 'p03'); subjectNames_recoded{idx} = 'Patient 3';
idx = contains(subjectNames, 'p04'); subjectNames_recoded{idx} = 'Patient 4';
idx = contains(subjectNames, 'p05'); subjectNames_recoded{idx} = 'Patient 5';
idx = contains(subjectNames, 'p06'); subjectNames_recoded{idx} = 'Patient 6';
idx = contains(subjectNames, 'p07'); subjectNames_recoded{idx} = 'Patient 7';
idx = contains(subjectNames, 'p10'); subjectNames_recoded{idx} = 'Patient 10';
idx = contains(subjectNames, 'p11'); subjectNames_recoded{idx} = 'Patient 11';

cmap = brewermap(length(subjectNames),'Set2');
%cmap(:,[2 3]) = 0.5;

nCols = 3;
ROW_IDX = [1 3 2]; % determines order of parameters

for jj = 1:nCols

	subplot(1,nCols,jj);hold on

    for ii = 1:length(subjectNames)
        
        idx = contains(D.channels.subject_name, subjectNames{ii});
        [~, ~, group_prob] = groupElecsByVisualArea(D.channels(idx,:), 'probabilisticresample');  
        [m_dat, se_dat] = averageWithinArea(derivedPrm(:,idx), group_prob, fun, nboot);
        %[m_mod, se_mod] = averageWithinArea(results.derived.params(:,idx), group_prob, fun, nboot);
        
        ix = find(~isnan(m_dat(1,:)));
        [hp, hc] = tde_plotPoints(m_dat(ROW_IDX(jj),ix)', squeeze(se_dat(ROW_IDX(jj),ix,:)), x(ix), 'errbar', 0, [], 20);
        hc.FaceColor = cmap(ii,:);
        hp.Color = cmap(ii,:);
        xlim([0 max(x)+1]); 
        if ROW_IDX(jj) == 1, ylim([0 0.5]), end
        if ROW_IDX(jj) == 2, ylim([0 1]), end
        if ROW_IDX(jj) == 3, ylim([0 1]), end

        ylabel(ylabels{ROW_IDX(jj)})

       
        set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
        xlabel('Visual area');
       
    end

    [~, channels, group_prob] = groupElecsByVisualArea(D.channels, 'probabilisticresample');   
    [m, se] = averageWithinArea(derivedPrm(ROW_IDX(jj),:), group_prob, [], 1000);
    h = tde_plotPoints(m', se, x, 'errbar', 0, '-');
    h.LineWidth = 2;
    h.MarkerSize = 30;
    if jj == 1
        legend([subjectNames_recoded; 'all']); legend box off
    end
  
end
set(findall(gcf,'-property','FontSize'),'FontSize',24)
