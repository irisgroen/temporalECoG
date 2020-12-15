channelsPRF = tde_getPRFparams(channels2fit);

figure;histogram(channelsPRF.aprf_ecc(chan_idx_R2),100)

% define cut offs
R2thresh    = 30;
eccthresh   = 3;
eccmax      = 10;

% select channels
chan_idx_R2     = channelsPRF.aprf_R2 > R2thresh;
chan_idx_eccfov = channelsPRF.aprf_ecc < eccthresh;
chan_idx_eccper = channelsPRF.aprf_ecc >= eccthresh & channelsPRF.aprf_ecc < eccmax;
chan_idx_area   = matchAreaNameToAtlas({'V123'}, channelsPRF.benson14_varea);

% combine criteria
fov_idx = chan_idx_R2 & chan_idx_eccfov & chan_idx_area;
per_idx = chan_idx_R2 & chan_idx_eccper & chan_idx_area;

% get the data
tmpdata_fov = data2fit(:,:,fov_idx);
tmpdata_per = data2fit(:,:,per_idx);

% plot
figure;
subplot(2,1,1);hold on
plot(flatten(mean(tmpdata_fov,3)),'k', 'LineWidth',2); 
plot(flatten(mean(tmpdata_per,3)),'r', 'LineWidth',2); 

set(gca, 'XTick',1:size(data2fit,1):size(data2fit,2)*size(data2fit,1), 'XTickLabel', []);
axis tight

l1 = sprintf('foveal (< %d degrees, n = %d)', eccthresh, length(find(fov_idx)));
l2 = sprintf('peripheral (> %d degrees, n = %d)', eccthresh, length(find(per_idx)));
legend({l1,l2}, 'Location', 'NorthWest');
title(sprintf('V123 (PRF R2 > %d)', R2thresh));


chan_idx_area   = matchAreaNameToAtlas({'higher'}, channelsPRF.benson14_varea);
fov_idx = chan_idx_R2 & chan_idx_eccfov & chan_idx_area;
per_idx = chan_idx_R2 & chan_idx_eccper & chan_idx_area;

tmpdata_fov = data2fit(:,:,fov_idx);
tmpdata_per = data2fit(:,:,per_idx);

subplot(2,1,2); hold on
plot(flatten(mean(tmpdata_fov,3)),'k', 'LineWidth',2); 
plot(flatten(mean(tmpdata_per,3)),'r', 'LineWidth',2); 

set(gca, 'XTick',1:size(data2fit,1):size(data2fit,2)*size(data2fit,1), 'XTickLabel', []);
axis tight
l1 = sprintf('foveal (< %d degrees, n = %d)', eccthresh, length(find(fov_idx)));
l2 = sprintf('peripheral (> %d degrees, n = %d)', eccthresh, length(find(per_idx)));
legend({l1,l2}, 'Location', 'NorthWest');
title(sprintf('higher (PRF R2 > %d)', R2thresh));

set(gcf, 'Position', get(0,'screensize'));
%set(gcf, 'Position', [1         436        1440         369]);