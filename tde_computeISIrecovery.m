function [ISIrecover] = tde_computeISIrecovery(data,t,stim_info)


[nSamp,~,nDatasets] = size(data);

% Set parameters

w = 0.3; % window for computing sum
shift = 0;
stimdur = 0.133;
conditionsOfInterest = {'ONEPULSE-4','ONEPULSE-5','TWOPULSE'};

ISIrecover = [];

for kk = 1:nDatasets

    % Get the data
    stim_idx = find(contains(stim_info.name, conditionsOfInterest));
    D = data(:,stim_idx,kk);
    stimNames = stim_info.name(stim_idx)';

    % Get onsets for each pulse for each stimulus
    nStim = length(stim_idx);
    pulse1_onset = zeros(nStim,1);
    pulse2_onset = stim_info.ISI(stim_idx) + stimdur;
    pulse2_onset(1) = w; % ONEPULSE-4 has no second stimulus

    % Compute response to first stimulus
    pulse1 = D;
    for ii = 1:nStim
        t_idx = t <= pulse2_onset(ii) + shift;
        pulse1(~t_idx,ii) = nan;
    end

    pulse1_mn = mean(pulse1,2,'omitnan');
    %pulse1_mn = pulse1(:,1); 
    %pulse1_mn = mean(pulse1(:,2:end),2,'omitnan');
    %pulse1_mn = mean(pulse1(:,3:end),2,'omitnan');
    %pulse1_mn = mean(pulse1(:,[end]),2,'omitnan');
    
    % figure;hold on;
    % plot(t,pulse1, 'LineWidth', 1);
    % plot(t,pulse1_mn, 'k', 'LineWidth', 2);
    % legend([stimNames 'mean']);

    % Compute response to second stimulus
    pulse2 = D;

    % Omit ONEPULSE-4 which has no second stimulus
    pulse2 = pulse2(:,2:end); 
    stimNames = stimNames(2:end);
    pulse2_onset = pulse2_onset(2:end);
    nStim = length(stimNames);

    % Subtract the mean response to pulse 1 from pulse 2
    pulse1_mn_to_subtract = pulse1_mn;
    pulse1_mn_to_subtract(isnan(pulse1_mn)) = 0;
    pulse2_sub = pulse2 - pulse1_mn_to_subtract;

    % figure;hold on;
    % plot(t,pulse2_sub, 'LineWidth', 2);
    % legend(stimNames);

    % Compute the sum over mean of pulse 1
    t_idx1 = t > pulse1_onset(1) + shift & t<= pulse1_onset(1) + w + shift;
    %pulse1_mn_summed = sum(pulse1_mn(t_idx1), 'omitnan');
    pulse1_mn_summed = max(pulse1_mn(t_idx1));
    
    % Compute the sum over each second pulse 
    t_idx2 = [];
    for ii = 1:nStim
        t_idx2(:,ii) = t > pulse2_onset(ii) + shift & t<= pulse2_onset(ii) + w + shift;
    end
    t_idx2 = logical(t_idx2);

    % Sum pulse 2
    pulse2_to_sum = pulse2_sub;
    pulse2_to_sum(~t_idx2) = 0;
    %pulse2_summed = sum(pulse2_to_sum,1, 'omitnan');
    pulse2_summed = max(pulse2_to_sum,[],1, 'omitnan');

    % Compute recovery
    ISIrecover(:,kk) = (pulse2_summed./pulse1_mn_summed) * 100; % in percentage

%     % Debug:
% 
%     % Plot pulse 1 mean + window
%     figure('Position', get(0, 'ScreenSize'));hold on;
%     subplot(round(nStim/2),2,1); hold on;
%     tmp = zeros(nSamp,1);
%     tmp(t_idx1) = 10;
%     plot(t,pulse1_mn, 'r','LineWidth', 2); 
%     plot(t,tmp, 'k');
%     title('mean pulse 1');
% 
%      % Plot pulse 2 + window
%     for ii = 1:nStim
%         tmp = zeros(nSamp,1);
%         tmp(t_idx2(:,ii)) = 10;
%         subplot(round(nStim/2),2, ii+1); hold on
%         plot(t,pulse2(:,ii),'b', 'LineWidth', 2);
%         plot(t,pulse2_sub(:,ii),'m:', 'LineWidth', 2);
%         plot(t,tmp, 'k');    
%         plot(t,pulse1_mn, 'r', 'LineWidth', 2);
%         title(stimNames{ii});
%     end
%     
%    % Plot recovery
%    %figure;plot(ISIrecover, 'k.-', 'MarkerSize', 50, 'LineWidth', 2);
    
%     % Debug 2 
%     
%     figure;hold on
%     subplot(1,3,1);hold on;
%     colors = parula(size(pulse2_sub,2));
% 
%     for ii = 1:size(pulse2_sub,2)
%         plot(t,pulse2_to_sum(:,ii), 'Color', colors(ii,:), 'LineWidth', 2);
%     end
%     plot(t,pulse1_mn, 'k', 'LineWidth', 3)
% 
%     dummy = nan(length(t),1);
%     dummy(t_idx1) = -10;
%     hold on
%     plot(t,dummy, 'k', 'LineWidth', 3)
% 
%     for ii = 1:size(t_idx2,2)
%         dummy = nan(length(t),1);
%         dummy(t_idx2(:,ii)) = -10+ii;
%         plot(t,dummy, 'Color', colors(ii,:), 'LineWidth', 3)
%     end
% 
%     subplot(1,3,2);hold on;
%     plot(pulse1_mn(t_idx1), 'k', 'LineWidth', 3)
%     for ii = 1:size(pulse2_sub,2)
%         tmp = pulse2_to_sum(t_idx2(:,ii),ii);
%         plot(tmp, 'Color', colors(ii,:), 'LineWidth', 2);
%     end
%     
%     subplot(1,3,3);hold on;
%     x = 1:size(pulse2_sub,2);
%     plot(1,pulse1_mn_summed, 'k.', 'MarkerSize', 100, 'LineStyle', 'none');
%     for ii = 1:length(x)
%         plot(x(ii)+1,pulse2_summed(ii), 'Color', colors(ii,:), 'Marker', '.', 'MarkerSize', 100, 'LineStyle', 'none')
%     end
%     xlim([0 max(x)+1]);
%     
% 	
%     set(gcf, 'Position', [41         303        1301         502]);
     
end

end
