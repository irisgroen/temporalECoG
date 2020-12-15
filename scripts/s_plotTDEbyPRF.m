channelsPRF = tde_getPRFparams(channels2fit);
channelsPRF.benson14_angle = channelsPRF.benson14_angle+90;

% define cut offs
R2thresh         = 30;
eccfovthresh     = 1.5;
eccparafovthresh = 6;
eccmax           = 12;

% select channels
chan_idx_R2         = channelsPRF.aprf_R2 > R2thresh;
chan_idx_eccfov     = channelsPRF.aprf_ecc < eccfovthresh;
chan_idx_eccparafov = channelsPRF.aprf_ecc >= eccfovthresh & channelsPRF.aprf_ecc < eccparafovthresh;
chan_idx_eccper     = channelsPRF.aprf_ecc >= eccparafovthresh & channelsPRF.aprf_ecc < eccmax;


%% contrast temporal time courses

% V123 %%%
chan_idx_area       = matchAreaNameToAtlas({'V123'}, channelsPRF.benson14_varea);

% combine criteria
fov_idx      = chan_idx_R2 & chan_idx_eccfov & chan_idx_area;
parafov_idx  = chan_idx_R2 & chan_idx_eccparafov & chan_idx_area;
per_idx      = chan_idx_R2 & chan_idx_eccper & chan_idx_area;

% get the data
tmpdata_fov     = data2fit(:,:,fov_idx);
tmpdata_parafov = data2fit(:,:,parafov_idx);
tmpdata_per     = data2fit(:,:,per_idx);

% plot
figure; hold on
subplot(2,1,1);hold on
plot(flatten(median(tmpdata_fov,3)),'k', 'LineWidth',2); 
plot(flatten(median(tmpdata_parafov,3)),'b', 'LineWidth',2); 
plot(flatten(mean(tmpdata_per,3)),'r', 'LineWidth',2); 

set(gca, 'XTick',1:size(data2fit,1):size(data2fit,2)*size(data2fit,1), 'XTickLabel', []);
axis tight

l1 = sprintf('foveal (<%0.1f degrees, n = %d)', eccfovthresh, length(find(fov_idx)));
l2 = sprintf('parafoveal (>%0.1f degrees, <%d degrees, n = %d)', eccfovthresh, eccparafovthresh, length(find(parafov_idx)));
l3 = sprintf('peripheral (>%d degrees, <%d degrees, n = %d)', eccparafovthresh, eccmax, length(find(per_idx)));
legend({l1,l2,l3}, 'Location', 'NorthWest');
title(sprintf('V123 (PRF R2 > %d)', R2thresh));

%%% HIGHER %%%
chan_idx_area   = matchAreaNameToAtlas({'higher'}, channelsPRF.benson14_varea);

% combine criteria
fov_idx      = chan_idx_R2 & chan_idx_eccfov & chan_idx_area;
parafov_idx  = chan_idx_R2 & chan_idx_eccparafov & chan_idx_area;%&~contains(channelsPRF.subject_name, 'som708');
per_idx      = chan_idx_R2 & chan_idx_eccper & chan_idx_area;

% get the data
tmpdata_fov     = data2fit(:,:,fov_idx);
tmpdata_parafov = data2fit(:,:,parafov_idx);
tmpdata_per     = data2fit(:,:,per_idx);

% plot
subplot(2,1,2);hold on
plot(flatten(median(tmpdata_fov,3)),'k', 'LineWidth',2); 
plot(flatten(median(tmpdata_parafov,3)),'b', 'LineWidth',2); 
plot(flatten(mean(tmpdata_per,3)),'r', 'LineWidth',2); 

set(gca, 'XTick',1:size(data2fit,1):size(data2fit,2)*size(data2fit,1), 'XTickLabel', []);
axis tight

l1 = sprintf('foveal (<%0.1f degrees, n = %d)', eccfovthresh, length(find(fov_idx)));
l2 = sprintf('parafoveal (>%0.1f degrees, <%d degrees, n = %d)', eccfovthresh, eccparafovthresh, length(find(parafov_idx)));
l3 = sprintf('peripheral (>%d degrees, <%d degrees, n = %d)', eccparafovthresh, eccmax, length(find(per_idx)));
legend({l1,l2,l3}, 'Location', 'NorthWest');
title(sprintf('higher (PRF R2 > %d)', R2thresh));

set(gcf, 'Position', get(0,'screensize'));
%set(gcf, 'Position', [1         436        1440         369]);

%% model params

objFunction = @LINEAR_RECTF_EXP_NORM_DELAY;
tmp = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', func2str(objFunction))));
paramNames = strsplit(tmp.params, ',');

% V123 %%%
chan_idx_area       = matchAreaNameToAtlas({'V123'}, channelsPRF.benson14_varea);

% combine criteria
fov_idx      = chan_idx_R2 & chan_idx_eccfov & chan_idx_area;
parafov_idx  = chan_idx_R2 & chan_idx_eccparafov & chan_idx_area;
per_idx      = chan_idx_R2 & chan_idx_eccper & chan_idx_area;

% get the data
tmpdata_fov     = params(:,fov_idx);
tmpdata_parafov = params(:,parafov_idx);
tmpdata_per     = params(:,per_idx);

% plot
figure; hold on
for p = 1:length(paramNames)
    subplot(2,length(paramNames),p);hold on;
    %plot(1, median(tmpdata_fov(p,:),2),'.k','MarkerSize', 30); 
    errorbar(0.95, mean(tmpdata_fov(p,:),2), std(tmpdata_fov(p,:),0,2)/sqrt(size(tmpdata_fov,2)), '.k', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
    %plot(1.05, median(tmpdata_parafov(p,:),2),'.b','MarkerSize', 30); 
	errorbar(1, mean(tmpdata_parafov(p,:),2), std(tmpdata_parafov(p,:),0,2)/sqrt(size(tmpdata_parafov,2)), '.b', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
    %plot(0.95, median(tmpdata_per(p,:),2),'.r', 'MarkerSize', 30); 
	%errorbar(1.05, median(tmpdata_per(p,:),2), std(tmpdata_per(p,:),0,2)/sqrt(size(tmpdata_per,2)), '.r', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
    title(paramNames{p});
    set(gca, 'XLim',[0.8 1.2], 'XTickLabel', []);
end

% V123 %%%
chan_idx_area       = matchAreaNameToAtlas({'higher'}, channelsPRF.benson14_varea);

% combine criteria
fov_idx      = chan_idx_R2 & chan_idx_eccfov & chan_idx_area;
parafov_idx  = chan_idx_R2 & chan_idx_eccparafov & chan_idx_area;
per_idx      = chan_idx_R2 & chan_idx_eccper & chan_idx_area;

% get the data
tmpdata_fov     = params(:,fov_idx);
tmpdata_parafov = params(:,parafov_idx);
tmpdata_per     = params(:,per_idx);

% plot
for p = 1:length(paramNames)
    subplot(2,length(paramNames),p+length(paramNames));hold on;
    %plot(1, median(tmpdata_fov(p,:),2),'.k','MarkerSize', 30); 
    errorbar(0.95, mean(tmpdata_fov(p,:),2), std(tmpdata_fov(p,:),0,2)/sqrt(size(tmpdata_fov,2)), '.k', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
    %plot(1.05, median(tmpdata_parafov(p,:),2),'.b','MarkerSize', 30); 
	errorbar(1,mean(tmpdata_parafov(p,:),2), std(tmpdata_parafov(p,:),0,2)/sqrt(size(tmpdata_parafov,2)), '.b', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
    %plot(0.95, median(tmpdata_per(p,:),2),'.r', 'MarkerSize', 30); 
	%errorbar(1.05, median(tmpdata_per(p,:),2), std(tmpdata_per(p,:),0,2)/sqrt(size(tmpdata_per,2)), '.r', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
    title(paramNames{p});
    set(gca, 'XLim',[0.8 1.2], 'XTickLabel', []);
end

set(gcf, 'Position', get(0,'screensize'));
%set(gcf, 'Position', [1         436        1440         369]);


%%
l1 = sprintf('foveal (<%0.1f degrees, n = %d)', eccfovthresh, length(find(fov_idx)));
l2 = sprintf('parafoveal (>%0.1f degrees, <%d degrees, n = %d)', eccfovthresh, eccparafovthresh, length(find(parafov_idx)));
l3 = sprintf('peripheral (>%d degrees, <%d degrees, n = %d)', eccparafovthresh, eccmax, length(find(per_idx)));
legend({l1,l2,l3}, 'Location', 'NorthWest');
title(sprintf('V123 (PRF R2 > %d)', R2thresh));

%%% HIGHER %%%
chan_idx_area   = matchAreaNameToAtlas({'higher'}, channelsPRF.benson14_varea);

% combine criteria
fov_idx      = chan_idx_R2 & chan_idx_eccfov & chan_idx_area;
parafov_idx  = chan_idx_R2 & chan_idx_eccparafov & chan_idx_area;%&~contains(channelsPRF.subject_name, 'som708');
per_idx      = chan_idx_R2 & chan_idx_eccper & chan_idx_area;

% get the data
tmpdata_fov     = data2fit(:,:,fov_idx);
tmpdata_parafov = data2fit(:,:,parafov_idx);
tmpdata_per     = data2fit(:,:,per_idx);

% plot
subplot(2,1,2);hold on
plot(flatten(median(tmpdata_fov,3)),'k', 'LineWidth',2); 
plot(flatten(median(tmpdata_parafov,3)),'b', 'LineWidth',2); 
plot(flatten(median(tmpdata_per,3)),'r', 'LineWidth',2); 

set(gca, 'XTick',1:size(data2fit,1):size(data2fit,2)*size(data2fit,1), 'XTickLabel', []);
axis tight

l1 = sprintf('foveal (<%0.1f degrees, n = %d)', eccfovthresh, length(find(fov_idx)));
l2 = sprintf('parafoveal (>%0.1f degrees, <%d degrees, n = %d)', eccfovthresh, eccparafovthresh, length(find(parafov_idx)));
l3 = sprintf('peripheral (>%d degrees, <%d degrees, n = %d)', eccparafovthresh, eccmax, length(find(per_idx)));
legend({l1,l2,l3}, 'Location', 'NorthWest');
title(sprintf('higher (PRF R2 > %d)', R2thresh));

set(gcf, 'Position', get(0,'screensize'));
%set(gcf, 'Position', [1         436        1440         369]);

%% comparison with benson

chan_idx_area       = matchAreaNameToAtlas({'V123'}, channelsPRF.benson14_varea);
chan_idx_ecc        = channelsPRF.aprf_ecc < 10;

idx                 = chan_idx_R2 & chan_idx_area & chan_idx_ecc;

figure;hold on;
subplot(2,3,1); hold on;
scatter(channelsPRF.benson14_angle(idx), channelsPRF.aprf_ang(idx), 100, 'k','filled');
title(sprintf('V123 angle: rho = %0.2f', corr(channelsPRF.benson14_angle(idx),channelsPRF.aprf_ang(idx), 'Type', 'Spearman', 'Rows', 'complete')));
xlabel('benson');ylabel('analyzePRF');
axis equal; set(gca, 'XLim', [0 360], 'YLim', [0 360]);
subplot(2,3,2); hold on;
scatter(channelsPRF.benson14_eccen(idx), channelsPRF.aprf_ecc(idx),  100, 'k','filled');
title(sprintf('V123 eccen: rho = %0.2f', corr(channelsPRF.benson14_eccen(idx),channelsPRF.aprf_ecc(idx), 'Type', 'Spearman', 'Rows', 'complete')));
axis equal; set(gca, 'XLim', [0 10], 'YLim', [0 10]);
subplot(2,3,3); hold on;
scatter(channelsPRF.benson14_sigma(idx), channelsPRF.aprf_rfsize(idx), 100, 'k', 'filled');
title(sprintf('V123 sigma: rho = %0.2f', corr(channelsPRF.benson14_sigma(idx),channelsPRF.aprf_rfsize(idx),'Type', 'Spearman',  'Rows', 'complete')));
axis equal; set(gca, 'XLim', [0 4], 'YLim', [0 4]);
%set(gcf, 'Position', [1         436        1440         369]);

chan_idx_area   = matchAreaNameToAtlas({'higher'}, channelsPRF.benson14_varea);

idx             = chan_idx_R2 & chan_idx_area & chan_idx_ecc;


subplot(2,3,4); hold on;
scatter(channelsPRF.benson14_angle(idx), channelsPRF.aprf_ang(idx), 100, 'k','filled');
title(sprintf('higher angle: rho = %0.2f', corr(channelsPRF.benson14_angle(idx),channelsPRF.aprf_ang(idx),'Type', 'Spearman',  'Rows', 'complete')));
xlabel('benson');ylabel('analyzePRF'); 
axis equal; set(gca, 'XLim', [0 360], 'YLim', [0 360]);
subplot(2,3,5); hold on;
scatter(channelsPRF.benson14_eccen(idx), channelsPRF.aprf_ecc(idx),  100, 'k','filled');
title(sprintf('higher eccen: rho = %0.2f', corr(channelsPRF.benson14_eccen(idx),channelsPRF.aprf_ecc(idx),'Type', 'Spearman',  'Rows', 'complete')));
axis equal; set(gca, 'XLim', [0 10], 'YLim', [0 10]);
subplot(2,3,6); hold on;
scatter(channelsPRF.benson14_sigma(idx), channelsPRF.aprf_rfsize(idx), 100, 'k', 'filled');
title(sprintf('higher sigma: rho = %0.2f', corr(channelsPRF.benson14_sigma(idx),channelsPRF.aprf_rfsize(idx), 'Type', 'Spearman', 'Rows', 'complete')));
axis equal; set(gca, 'XLim', [0 8], 'YLim', [0 8]);

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
