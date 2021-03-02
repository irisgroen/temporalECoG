
% individual electrodes and DN model fits
modelfun = @LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datatype = 'individualelecs';
[d] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

%%
channelsPRF = tde_getPRFparams(d.channels);
[channelsPRF] = convertBensonangleToAnalyzePRFangle(channelsPRF);

% define cut offs
R2thresh         = 30;
eccfovthresh     = 2;
eccparafovthresh = 16;
eccmax           = 20;

% select channels
chan_idx_R2         = channelsPRF.aprf_R2 > R2thresh;
chan_idx_eccfov     = channelsPRF.aprf_ecc < eccfovthresh;
chan_idx_eccparafov = channelsPRF.aprf_ecc >= eccfovthresh & channelsPRF.aprf_ecc < eccparafovthresh;
chan_idx_eccper     = channelsPRF.aprf_ecc >= eccparafovthresh & channelsPRF.aprf_ecc < eccmax;

%% contrast temporal time courses
%areaNames = {'V1','V2','V3','higher'};
areaNames = {'V123','higher'};
fun = @mean;
numboot = 1000;

[~, ~, group_prob] = groupElecsByVisualArea(d.channels, 'probabilisticresample', areaNames );   

% combine criteria
fov_idx      = chan_idx_R2 & chan_idx_eccfov; %&~contains(channelsPRF.subject_name, 'som708'); 
parafov_idx  = chan_idx_R2 & chan_idx_eccparafov; %&~contains(channelsPRF.subject_name, 'som708'); 
per_idx      = chan_idx_R2 & chan_idx_eccper;% &~contains(channelsPRF.subject_name, 'som708');

% average the electrodes
[tmpdata_fov] = averageWithinArea(d.data(:,:,fov_idx), group_prob(fov_idx,:), fun, numboot);
[tmpdata_parafov] = averageWithinArea(d.data(:,:,parafov_idx), group_prob(parafov_idx,:), fun, numboot);
[tmpdata_per] = averageWithinArea(d.data(:,:,per_idx), group_prob(per_idx,:), fun, numboot);

% plot
figure; hold on

for ii = 1:length(areaNames)
    subplot(length(areaNames),1,ii);hold on
    toplot = flatten(tmpdata_fov(:,:,ii));
    plot(smooth(toplot/norm(toplot,2),20),'k', 'LineWidth',2); 
    toplot = flatten(tmpdata_parafov(:,:,ii));
    plot(smooth(toplot/norm(toplot,2),20),'r', 'LineWidth',2); 
    %plot(smooth(flatten(tmpdata_per(:,:,ii)),20),'b', 'LineWidth',2); 

    set(gca, 'XTick',1:size(d.data,1):size(d.data,2)*size(d.data,1), 'XTickLabel', []);
    axis tight

    l1 = sprintf('foveal (<%0.1f degrees, n = %d)', eccfovthresh, length(find(group_prob(fov_idx,ii)>0)));
    l2 = sprintf('peripheral (>%0.1f degrees, <%d degrees, n = %d)', eccfovthresh, eccparafovthresh, length(find(group_prob(parafov_idx,ii)>0)));
    %l3 = sprintf('peripheral (>%d degrees, <%d degrees, n = %d)', eccparafovthresh, eccmax, length(find(group_prob(per_idx,ii)>0)));
    legend({l1,l2}, 'Location', 'NorthWest');
    %legend({l1,l2,l3}, 'Location', 'NorthWest');
    title(sprintf('%s (PRF R2 > %d)', areaNames{ii}, R2thresh));
end
set(gcf, 'Position', get(0,'screensize'));
%set(gcf, 'Position', [1         436        1440         369]);

%% model params

objFunction = @LINEAR_RECTF_EXP_NORM_DELAY;
tmp = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', func2str(objFunction))));
paramNames = strsplit(tmp.params, ',');

%areaNames = {'V1','V2','V3','higher'};
areaNames = {'V123','higher'};
fun = @median;
numboot = 10000;

[~, ~, group_prob] = groupElecsByVisualArea(d.channels, 'probabilisticresample', areaNames );   

% combine criteria
fov_idx      = chan_idx_R2 & chan_idx_eccfov;
parafov_idx  = chan_idx_R2 & chan_idx_eccparafov;
per_idx      = chan_idx_R2 & chan_idx_eccper;

% average the electrodes
[d_fov, se_fov] = averageWithinArea(d.params(:,fov_idx), group_prob(fov_idx,:), fun, numboot);
[d_pfov, se_pfov] = averageWithinArea(d.params(:,parafov_idx), group_prob(parafov_idx,:), fun, numboot);
%[d_per, se_per] = averageWithinArea(params(:,:,per_idx), group_prob(per_idx,:), fun, numboot);

% plot
figure; hold on
cp = 1;
for ii = 1:length(areaNames)
    for p = 1:length(paramNames)
        subplot(length(areaNames),length(paramNames),cp);hold on;
        errorbar(0.95, d_fov(p,ii), d_fov(p,ii)-se_fov(p,ii,1), se_fov(p,ii,2)-d_fov(p,ii), '.k', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
        errorbar(1, d_pfov(p,ii), d_pfov(p,ii)-se_pfov(p,ii,1), se_pfov(p,ii,2)-d_pfov(p,ii), '.r', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
        %errorbar(1.05, d_per(p,ii), d_per(p,ii)-se_per(p,ii,1), se_per(p,ii,2)-d_per(p,ii), '.b', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
        title(sprintf('%s %s', areaNames{ii}, paramNames{p}));
        set(gca, 'XLim',[0.8 1.2], 'XTickLabel', []);
        cp = cp + 1;
    end
end
set(gcf, 'Position', get(0,'screensize'));
%        set(gcf, 'Position', [1         436        1440         369]);


%% correlations
%corr(channelsPRF.aprf_ecc(chan_idx_R2 &  chan_idx_area), params(4,chan_idx_R2 &  chan_idx_area)', 'Type', 'Spearman')
areaName = {'higher'};
chan_idx_area       = matchAreaNameToAtlas(areaName, channelsPRF.benson14_varea);
subjectNames        = unique(channelsPRF.subject_name(chan_idx_R2 & chan_idx_area));
colors              = 1-parula(length(subjectNames));

figure;hold on
for p = 1:length(paramNames)
	subplot(2,5,p);hold on;
    for ii = 1:length(subjectNames)
        idx = chan_idx_R2 & contains(channelsPRF.subject_name, subjectNames{ii}) & chan_idx_area;
        scatter(channelsPRF.aprf_ecc(idx), d.params(p,idx), 100, colors(ii,:), 'filled')
    end
    xlabel('aPRF eccentricity');
    title(sprintf('%s %s', areaName{:}, paramNames{p}));
    set(gca, 'XLim', [0 10]);
end
legend(subjectNames);
set(gcf, 'Position', get(0,'screensize'));

%% comparison with benson

chan_idx_ecc        = channelsPRF.aprf_ecc < 10;
idx                 = chan_idx_R2 & chan_idx_ecc;

figure;hold on;
line([0 360], [0 360], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
scatter(channelsPRF.benson14_angle(idx), channelsPRF.aprf_ang(idx), 100, 'k','filled');
title(sprintf('V123 angle: circular r = %0.2f', circ_corrcc(circ_ang2rad(channelsPRF.benson14_angle(idx)),circ_ang2rad(channelsPRF.aprf_ang(idx)))));
xlabel('benson');ylabel('analyzePRF');

%%
chan_idx_area       = matchAreaNameToAtlas({'V123'}, channelsPRF.benson14_varea);
chan_idx_ecc        = channelsPRF.aprf_ecc < 20;

idx                 = chan_idx_R2 & chan_idx_area & chan_idx_ecc;

figure;hold on;
subplot(2,3,1); hold on;
line([0 360], [0 360], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
scatter(channelsPRF.benson14_angle(idx), channelsPRF.aprf_ang(idx), 100, 'k','filled');
title(sprintf('V123 angle: circular r = %0.2f', circ_corrcc(deg2rad(channelsPRF.benson14_angle(idx)),deg2rad(channelsPRF.aprf_ang(idx)))));
xlabel('benson');ylabel('analyzePRF');
axis equal; set(gca, 'XLim', [0 360], 'YLim', [0 360]);
subplot(2,3,2); hold on;
line([0 10], [0 10], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
scatter(channelsPRF.benson14_eccen(idx), channelsPRF.aprf_ecc(idx),  100, 'k','filled');
title(sprintf('V123 eccen: rho = %0.2f', corr(channelsPRF.benson14_eccen(idx),channelsPRF.aprf_ecc(idx), 'Type', 'Spearman', 'Rows', 'complete')));
axis equal; set(gca, 'XLim', [0 10], 'YLim', [0 10]);
subplot(2,3,3); hold on;
line([0 4], [0 4], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
scatter(channelsPRF.benson14_sigma(idx), channelsPRF.aprf_rfsize(idx), 100, 'k', 'filled');
title(sprintf('V123 sigma: rho = %0.2f', corr(channelsPRF.benson14_sigma(idx),channelsPRF.aprf_rfsize(idx),'Type', 'Spearman',  'Rows', 'complete')));
axis equal; set(gca, 'XLim', [0 3], 'YLim', [0 3]);
%set(gcf, 'Position', [1         436        1440         369]);

chan_idx_area   = matchAreaNameToAtlas({'higher'}, channelsPRF.benson14_varea);

idx             = chan_idx_R2 & chan_idx_area & chan_idx_ecc;


subplot(2,3,4); hold on;
line([0 360], [0 360], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
scatter(channelsPRF.benson14_angle(idx), channelsPRF.aprf_ang(idx), 100, 'k','filled');
title(sprintf('higher angle: circular r = %0.2f', circ_corrcc(deg2rad(channelsPRF.benson14_angle(idx)),deg2rad(channelsPRF.aprf_ang(idx)))));
xlabel('benson');ylabel('analyzePRF'); 
axis equal; set(gca, 'XLim', [0 360], 'YLim', [0 360]);
subplot(2,3,5); hold on;
line([0 10], [0 10], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
scatter(channelsPRF.benson14_eccen(idx), channelsPRF.aprf_ecc(idx),  100, 'k','filled');
title(sprintf('higher eccen: rho = %0.2f', corr(channelsPRF.benson14_eccen(idx),channelsPRF.aprf_ecc(idx),'Type', 'Spearman',  'Rows', 'complete')));
axis equal; set(gca, 'XLim', [0 10], 'YLim', [0 10]);
subplot(2,3,6); hold on;
line([0 8], [0 8], 'LineStyle', ':', 'Color', [0.5 0.5 0.5]);
scatter(channelsPRF.benson14_sigma(idx), channelsPRF.aprf_rfsize(idx), 100, 'k', 'filled');
title(sprintf('higher sigma: rho = %0.2f', corr(channelsPRF.benson14_sigma(idx),channelsPRF.aprf_rfsize(idx), 'Type', 'Spearman', 'Rows', 'complete')));
axis equal; set(gca, 'XLim', [0 6], 'YLim', [0 6]);

set(gcf, 'Position', get(0,'screensize'));

%% summary stats
%figure;histogram(channelsPRF.aprf_ecc(chan_idx_R2),100)

figure;hold on;

chan_idx_area   = matchAreaNameToAtlas({'V123'}, channelsPRF.benson14_varea);
idx             = chan_idx_R2 & chan_idx_area;

subplot(2,4,1); hold on;
histogram(channelsPRF.aprf_R2(idx),50);
xlabel('R2'); ylabel('# elecs');xlim([0 100]);title(sprintf('V123 (n = %d)', length(find(idx))));
subplot(2,4,2); hold on;
histogram(channelsPRF.aprf_ang(idx),50);
xlabel('angle'); xlim([0 360]);
subplot(2,4,3); hold on; xlim([0 25]);
histogram(channelsPRF.aprf_ecc(idx),50);
xlabel('eccentricity'); 
subplot(2,4,4); hold on; xlim([0 8]);
histogram(channelsPRF.aprf_rfsize(idx),50);
xlabel('rfsize');

chan_idx_area   = matchAreaNameToAtlas({'higher'}, channelsPRF.benson14_varea);
idx             = chan_idx_R2 & chan_idx_area;

subplot(2,4,5); hold on;
histogram(channelsPRF.aprf_R2(idx),50);
xlabel('R2'); ylabel('# elecs'); xlim([0 100]); title(sprintf('higher (n = %d)', length(find(idx))));
subplot(2,4,6); hold on;
histogram(channelsPRF.aprf_ang(idx),50);
xlabel('angle');  xlim([0 360]);
subplot(2,4,7); hold on;
histogram(channelsPRF.aprf_ecc(idx),50);
xlabel('eccentricity');  xlim([0 25]);
subplot(2,4,8); hold on;
histogram(channelsPRF.aprf_rfsize(idx),50);
xlabel('rfsize');xlim([0 8]);

set(gcf, 'Position', get(0,'screensize'));
