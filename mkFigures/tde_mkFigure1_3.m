% Figure 1_3 : 
% Fitted parameters for DN model
%
% 2022 Iris Groen
modelfun = @DN;
xvalmode = 1;
datatype = 'individualelecs';
[D] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

[results] = tde_evaluateModelFit(D,0);
numboot = 10000;

[~, channels, group_prob] = groupElecsByVisualArea(D.channels, 'probabilisticresample');
nChans = height(channels);

%% Extended data: fitted parameters

figure(1); clf; hold on
set(gcf, 'position',  get(0, 'screensize'));

% Read in parameter names from json
json = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', 'DN')));
paramNames = strsplit(json.params,',');
nParams = length(paramNames);

% Plot fitted parameters
for p = 1:nParams
    subplot(2,ceil(nParams/2),p); hold on    
    [m, se] = averageWithinArea(results.params(p,:), group_prob);
    errorbar(1:nChans, m, m-se(:,1)', se(:,2)'-m, '.k', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
    set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
    title(paramNames{p}); xlabel('Visual area'); ylabel('Parameter value');set(gca, 'fontsize', 16);
    if p == 1; ylim([0 0.1]); end
    if p == 2; ylim([0 0.6]); end
    if p == 3; ylim([0 0.4]); end
    if p == 4; ylim([0.8 1.8]); end
    if p == 5; ylim([0 0.1]); end
    if p == 6; ylim([0 0.15]); end
    if p == 7; ylim([0 4]); end
end

set(findall(gcf,'-property','FontSize'),'FontSize',20)


