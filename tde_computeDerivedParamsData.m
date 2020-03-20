function [params] = tde_computeDerivedParamsData(data,channels,t,shift,stim_info)

%load('/Volumes/server/Projects/BAIR/Papers/TemporalDynamicsECoG/results/DN_xvalmode0_20200312T134253.mat')

[~,~,nDatasets] = size(data);

stim_idx = find(contains(stim_info.name, {'ONEPULSE-4', 'TWOPULSE'}));
nStim = length(stim_idx);
stimdur = unique(stim_info.duration(stim_idx));
stimISI = stim_info.ISI(stim_idx);
%shift = ones(nDatasets,1)*0.05;

data_pulse1 = data(:,stim_idx,:);
data_pulse2 = data(:,stim_idx,:);

for ii = 1:nDatasets
     for ss = 1:nStim
         if ss == 1
             isi = stimISI(end);
         else
             isi = stimISI(ss);
         end
%          if ss > 4
%              isi = stimISI(4);
%          else
%              isi = stimISI(ss);
%          end
        pulse1_tidx = t > shift(ii) & t < shift(ii) + stimdur + isi;
        pulse2_tidx = t > shift(ii) + stimdur + isi;        
        data_pulse1(~pulse1_tidx,ss,ii) = nan;
        data_pulse2(~pulse2_tidx,ss,ii) = nan;

    end 
end

data_pulse1_mn = mean(data_pulse1,2,'omitnan');
%data_pulse1_mn(t>0.3)=nan;
data_pulse1_mn_to_subtract = data_pulse1_mn;
data_pulse1_mn_to_subtract(isnan(data_pulse1_mn)) = 0;
data_pulse2_sub = data_pulse2 - data_pulse1_mn_to_subtract;

pulse1_mn_summed = [];
pulse1_summed = [];
pulse2_summed = [];

% plots
cmap = num2cell(copper(length(stimISI)), 2);
for ii = 1:nDatasets

    figure('Name', channels.name{ii}, 'Position',  [78 61 1098 637]);
    subplot(2,3,1);
    h = plot(t,data(:,stim_idx,ii), 'LineWidth',2);
    set(h, {'color'}, cmap);
    ylims = get(gca, 'YLim'); set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
    title('selected conditions');
    xlabel('Time (s)');

    subplot(2,3,2);
    h = plot(t,data_pulse1(:,:,ii), 'LineWidth',2);
    set(h, {'color'}, cmap);
    set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
    xlim([t(1) t(end)]);
    title('response to first pulse');
    xlabel('Time (s)');
    
    subplot(2,3,3);
    h = plot(t,data_pulse1_mn(:,:,ii), 'k', 'LineWidth',2);
    set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
    xlim([t(1) t(end)]);
    title('response to first pulse averaged');
    xlabel('Time (s)');

    subplot(2,3,4);   
    h = plot(t,data_pulse2(:,:,ii), 'LineWidth',2);
    set(h, {'color'}, cmap);
    set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
    xlim([t(1) t(end)]);
    title('response to second pulse');
    xlabel('Time (s)');
    
    subplot(2,3,5);   
    h = plot(t,data_pulse2_sub(:,:,ii), 'LineWidth',2);
    set(h, {'color'}, cmap);
    set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
    xlim([t(1) t(end)]);
    title('response to second pulse - first pulse');
    xlabel('Time (s)');

    subplot(2,3,6);hold on
    pulse1_mn_summed(ii) = squeeze(sum(data_pulse1_mn(:,:,ii),1, 'omitnan'));
    pulse1_summed(ii) = squeeze(sum(data_pulse1(:,1,ii),1, 'omitnan'));
    pulse2_summed(:,ii) = squeeze(sum(data_pulse2_sub(:,:,ii),1, 'omitnan'));
    line([stimISI(1) stimISI(end)], [pulse1_mn_summed(ii) pulse1_mn_summed(ii)],'Color','r', 'LineStyle', '-', 'LineWidth',2)
	line([stimISI(1) stimISI(end)], [pulse1_summed(ii) pulse1_summed(ii)],'Color','r', 'LineStyle', '--', 'LineWidth',2)
    plot(stimISI,pulse2_summed(:,ii), '-ko','LineWidth',2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    title('summed response to second pulse');
    xlabel('ISI (s)');
    xlim([stimISI(1)-0.1 stimISI(end)+0.1]);
    legend({'first pulse mean', 'first pulse ONEPULSE-4', 'second pulse'}, 'Location', 'SouthEast')

end

figure; hold on;
h1 = scatter(ones(nDatasets,1)*-0.02,pulse1_summed, 150, parula(nDatasets),'filled');
set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
h = plot(stimISI,pulse2_summed, '-o','LineWidth',2, 'MarkerSize', 10);
xlabel('ISI (s)');
ylabel('summed response');
set(h, {'color'}, num2cell(parula(nDatasets), 2));
xlim([stimISI(1)-0.1 stimISI(end)+0.1]);
legend(channels.name);
title('recovery with ISI - all areas');

end