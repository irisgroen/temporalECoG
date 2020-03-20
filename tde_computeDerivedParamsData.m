function [params] = tde_computeDerivedParamsData(data,channels,t,shift,stim_info)

%load('/Volumes/server/Projects/BAIR/Papers/TemporalDynamicsECoG/results/DN_xvalmode0_20200312T134253.mat')

[~,~,nDatasets] = size(data);

stim_idx = find(contains(stim_info.name, {'ONEPULSE-4', 'TWOPULSE'}));
nStim = length(stim_idx);
stimdur = unique(stim_info.duration(stim_idx));
stimISI = stim_info.ISI(stim_idx);

%shift = ones(nDatasets,1)*0.05;

pulse1_summed = [];
pulse2_summed = [];


for ii = 1:nDatasets
    
    pulse1_tidx = t > shift(ii) & t < shift(ii) + stimdur;
    pulse2_tidx = t > shift(ii) + stimdur;

    data_tmp1 = data;
    data_tmp1(~pulse1_tidx,:,:) = nan;
    data_tmp2 = data;
    data_tmp2(~pulse2_tidx,:,:) = nan;
    
   
end

% plots
for ii = 1:nDatasets

    figure('Name', channels.name{ii}, 'Position',  [360    61   816   637]);
    subplot(2,2,1);
    h = plot(t,data(:,stim_idx,ii), 'LineWidth',2);
    set(h, {'color'}, num2cell(copper(length(stimISI)), 2));
    ylims = get(gca, 'YLim'); set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
    title('selected conditions');
    xlabel('Time (s)');

    subplot(2,2,2);
    h = plot(t,data_tmp1(:,stim_idx,ii), 'LineWidth',2);
    set(h, {'color'}, num2cell(copper(length(stimISI)), 2));
    set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
    xlim([t(1) t(end)]);
    title('response to first pulse');
    xlabel('Time (s)');

    subplot(2,2,3);   
    h = plot(t,data_tmp2(:,stim_idx,ii), 'LineWidth',2);
    set(h, {'color'}, num2cell(copper(length(stimISI)), 2));
    set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
    xlim([t(1) t(end)]);
    title('response to second pulse');
    xlabel('Time (s)');

    subplot(2,2,4);hold on
    pulse1_summed(ii) = mean(squeeze(sum(data_tmp1(:,stim_idx,ii),1, 'omitnan')),2);
    pulse2_summed(:,ii) = squeeze(sum(data_tmp2(:,stim_idx,ii),1, 'omitnan'));
    line([stimISI(1) stimISI(end)], [pulse1_summed(ii) pulse1_summed(ii)],'Color','r', 'LineStyle', '--', 'LineWidth',2)
    plot(stimISI,pulse2_summed(:,ii), '-ko','LineWidth',2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    title('summed response to second pulse');
    xlabel('ISI');
    xlim([stimISI(1)-0.1 stimISI(end)+0.1]);
    legend({'first pulse', 'second pulse'}, 'Location', 'SouthEast')

end

figure; hold on;
h1 = scatter(ones(nDatasets,1)*-0.02,pulse1_summed, 150, parula(nDatasets),'filled');
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
h = plot(stimISI,pulse2_summed, '-o','LineWidth',2, 'MarkerSize', 10);
xlabel('ISI');
ylabel('summed response');
set(h, {'color'}, num2cell(parula(nDatasets), 2));
xlim([stimISI(1)-0.1 stimISI(end)+0.1]);
legend(channels.name);
title('recovery with ISI - all areas');

end