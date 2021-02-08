function [ISIrecover,ts,w] = tde_computeISIrecovery_sim(data,t,stim_info,srate,w,shift,metric)
% Computes recovery without using the ONE PULSE conditions for arbitrary
% TWOPULSE names

if ~exist('w', 'var') || isempty(w)
    w = 0.3; % window for computing recover
end

if ~exist('shift', 'var') || isempty(shift)
    shift = 0; % shift of window relative to stimulus onset
end

if ~exist('metric', 'var') || isempty(metric)
    metric = 'sum'; % 'sum' or 'max'
end

[~,~,nDatasets] = size(data);

% Set parameters
conditionsOfInterest = {'TWOPULSE'};
stim_idx = find(contains(stim_info.name, conditionsOfInterest));
stimdur = stim_info.duration(stim_idx(1));
% Make sure the window size matches the sample rate
w = round(w*srate)*(1/srate);

ISIrecover = [];
ts = [];

for kk = 1:nDatasets

    % Get the data
    D = data(:,stim_idx,kk);

    % Get onsets for each pulse for each stimulus
    nStim = length(stim_idx);
    pulse1_onset = zeros(nStim,1);
    pulse2_onset = stim_info.ISI(stim_idx) + stimdur;
    
    % Compute response to first stimulus
    pulse1 = D;
    for ii = 1:nStim
        t_idx = t <= pulse2_onset(ii) + shift;
        pulse1(~t_idx,ii) = nan;
    end

    pulse1_mn = mean(pulse1,2,'omitnan');
    % figure;hold on;
    % plot(t,pulse1, 'LineWidth', 1);
    % plot(t,pulse1_mn, 'k', 'LineWidth', 2);
    % legend([stimNames 'mean']);

    % Compute response to second stimulus
    pulse2 = D;

    % Subtract the mean response to pulse 1 from pulse 2
    pulse1_mn_to_subtract = pulse1_mn;
    pulse1_mn_to_subtract(isnan(pulse1_mn)) = 0;
    pulse2_sub = pulse2 - pulse1_mn_to_subtract;

    % figure;hold on;
    % plot(t,pulse2_sub, 'LineWidth', 2);
    % legend(stimNames);

    % Compute the sum over mean of pulse 1
    t_idx1 = t > pulse1_onset(1) + shift & t<= pulse1_onset(1) + w + shift;
    switch metric
        case 'sum'
            pulse1_mn_summed = sum(pulse1_mn(t_idx1), 'omitnan');
        case 'max'
            pulse1_mn_summed = max(pulse1_mn(t_idx1));
        otherwise
            error('unknown metric, choose sum or max')
    end
    
    % Compute the sum over each second pulse 
    t_idx2 = [];
    for ii = 1:nStim
        t_idx2(:,ii) = t > pulse2_onset(ii) + shift & t<= pulse2_onset(ii) + w + shift;
    end
    t_idx2 = logical(t_idx2);

    % Sum pulse 2
    pulse2_to_sum = pulse2_sub;
    pulse2_to_sum(~t_idx2) = 0;
    switch metric
        case 'sum'
            pulse2_summed = sum(pulse2_to_sum,1, 'omitnan');
        case 'max'
            pulse2_summed = max(pulse2_to_sum,[],1, 'omitnan');
    end

    % Compute recovery
    ISIrecover(:,kk) = (pulse2_summed./pulse1_mn_summed); % in percentage
    ts(:,1,kk) = pulse1_mn(t_idx1);
    for ii = 1:size(pulse2_sub,2)
        ts(:,ii+1,kk) = pulse2_to_sum(t_idx2(:,ii),ii);
    end
    
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
%     % Plot recovery
%     figure;plot(ISIrecover, 'k.-', 'MarkerSize', 50, 'LineWidth', 2);
     
end

end
