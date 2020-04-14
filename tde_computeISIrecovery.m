function [ISIrecover] = tde_computeISIrecovery(data,channels,t,stim_info)



[~,~,nDatasets] = size(data);

% Trecover
stim_idx = find(contains(stim_info.name, {'ONEPULSE-4','ONEPULSE-5','TWOPULSE'}));
nStim = length(stim_idx);
stimdur = unique(stim_info.duration(stim_idx(1)));
stimISI = stim_info.ISI(stim_idx);
shift = ones(nDatasets,1)*0.05;

data_pulse1 = data(:,stim_idx,:);
data_pulse2 = data(:,stim_idx,:);

for ii = 1:nDatasets
     for ss = 1:nStim
         if ss == 1
             isi = stimISI(5);
         else
             isi = stimISI(ss);
         end
         if ss < 5
            pulse1_tidx = t > 0 & t < shift(ii) + stimdur + isi;
         else
            pulse1_tidx = t > 0 & t < shift(ii) + stimdur + isi;
         end
        pulse2_tidx = t > shift(ii) + stimdur + isi & t < shift(ii) + stimdur + isi + stimdur + 0.3;% shift(ii);        
        data_pulse1(~pulse1_tidx,ss,ii) = nan;
        data_pulse2(~pulse2_tidx,ss,ii) = nan;
    end 
end

% subtract first pulse from second pulse
data_pulse1_mn = mean(data_pulse1,2,'omitnan');
data_pulse1_mn_to_subtract = data_pulse1_mn;
data_pulse1_mn_to_subtract(isnan(data_pulse1_mn)) = 0;
data_pulse2_sub = data_pulse2 - data_pulse1_mn_to_subtract;
%data_pulse2_sub = data(:,stim_idx,:) - data_pulse1_mn_to_subtract;

% compute one pulse baseline
pulse1_mn_tidx = t>0&t<0.55;
pulse1_mn_summed = squeeze(sum(data_pulse1_mn(pulse1_mn_tidx,:,:),1, 'omitnan'))';
data_pulse1_mn(~pulse1_mn_tidx,:,:)=nan; % for plotting

% plots
cmap = num2cell(parula(length(stimISI)), 2);
pulse2_summed = squeeze(sum(data_pulse2_sub(:,2:end,:),1, 'omitnan'));

ISIrecover = pulse2_summed(2:end,:)./pulse1_mn_summed;

end
%pulse2_summed = [];
% for ii = 1:nDatasets
% 
%     figure('Name', channels.name{ii}, 'Position',  [78 61 1098 637]);
%     subplot(2,3,1);
%     h = plot(t,data(:,stim_idx,ii), 'LineWidth',2);
%     set(h, {'color'}, cmap);
%     ylims = get(gca, 'YLim'); set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
%     title('selected conditions');
%     xlabel('Time (s)');
% 
%     subplot(2,3,2);
%     h = plot(t,data_pulse1(:,:,ii), 'LineWidth',2);
%     set(h, {'color'}, cmap);
%     set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
%     xlim([t(1) t(end)]);
%     title('response to first pulse');
%     xlabel('Time (s)');
%     
%     subplot(2,3,3);
%     h = plot(t,data_pulse1_mn(:,:,ii), 'k', 'LineWidth',2);
%     set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
%     xlim([t(1) t(end)]);
%     title('response to first pulse averaged');
%     xlabel('Time (s)');
% 
%     subplot(2,3,4);   
%     h = plot(t,data_pulse2(:,2:end,ii), 'LineWidth',2);
%     set(h, {'color'}, cmap(2:end,:));
%     set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
%     xlim([t(1) t(end)]);
%     title('response to second pulse');
%     xlabel('Time (s)');
%     legend(stim_info.name(stim_idx(2:end)));
%     
%     subplot(2,3,5);   
%     h = plot(t,data_pulse2_sub(:,2:end,ii), 'LineWidth',2);
%     set(h, {'color'}, cmap(2:end,:));
%     set(gca, 'YLim',[-0.5 ylims(2)+0.1*ylims(2)]);
%     xlim([t(1) t(end)]);
%     title('response to second pulse - first pulse');
%     xlabel('Time (s)');
% 
%     subplot(2,3,6);hold on
%     %pulse1_summed(ii) = squeeze(sum(data_pulse1(:,1,ii),1, 'omitnan'));
%     pulse2_summed(:,ii) = squeeze(sum(data_pulse2_sub(:,2:end,ii),1, 'omitnan'));
%     %line([stimISI(1) stimISI(end)], [pulse1_mn_summed(ii) pulse1_mn_summed(ii)],'Color','r', 'LineStyle', '-', 'LineWidth',2)
% 	%line([stimISI(1) stimISI(end)], [pulse1_summed(ii) pulse1_summed(ii)],'Color','r', 'LineStyle', '--', 'LineWidth',2)
%     plot(stimISI(2:end),pulse2_summed(:,ii)/pulse1_mn_summed(ii), '-ko','LineWidth',2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
%     title('summed response to second pulse');
%     xlabel('ISI (s)');
%     xlim([stimISI(1)-0.1 stimISI(end)+0.1]);
%     ylim([0 1]);
%     %legend({'first pulse mean', 'first pulse ONEPULSE-4', 'second pulse'}, 'Location', 'SouthEast')
% 
% end
% 
% figure; hold on;
% %h1 = scatter(ones(nDatasets,1)*-0.02,pulse1_summed, 150, parula(nDatasets),'filled');
% %set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
% ind = nDatasets:-1:1;
% h = plot(stimISI(2:end),pulse2_summed(:,ind)./pulse1_mn_summed(ind), '.-','LineWidth',2, 'MarkerSize', 50);
% xlabel('ISI (s)');
% ylabel('summed response');
% set(h, {'color'}, num2cell(flipud(parula(nDatasets)), 2));
% xlim([stimISI(1)-0.1 stimISI(end)+0.1]);
% ylim([0 1]);
% legend(channels.name(ind), 'Location', 'SouthEast');
% title('recovery with ISI - all areas');

%end