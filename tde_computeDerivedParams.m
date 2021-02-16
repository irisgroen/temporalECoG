function [derivedPrm, paramNames, pred_sustained, pred_transient, t] = tde_computeDerivedParams(objFunction, prm)
% Simulate model response to a long sustained stimulus
% Extract four summary statistics:
% - t2pk    = time to peak 
% - r_asymp = magnitude or response at asymptote
% - Rdouble = changed in summed response with doubling of duration
% - t_isi = magnitude or response at asymptote

paramNames = {'t2pk', 'rAsymp', 'wSize'};
derivedPrm = [];

%% Compute derived parameters Time2Peak and Rasymptote

% Simulate response to a long stimulus
t    = 0.001 : 0.001 : 2;
stim = ones(length(t),1);
srate = 1/median(diff(t)); % samples per second

[~, pred] = objFunction(prm, [], stim, srate);

[~,x] = max(pred);
% Find the timepoint corresponding to the maximum magnitude
derivedPrm(1)    = t(x);
% Find the ratio between the asymptote and the maximum magnitude
derivedPrm(2)    = pred(end)/max(pred);
pred_sustained   = pred;

%% Compute derived parameters WindowSize

% Simulate response to a short (5 ms) pulse 
t    = 0.001 : 0.001 : 0.5;
stim = zeros(length(t),1);
stim(1:5) = 1;
srate = 1/median(diff(t)); % samples per second

[~, pred] = objFunction(prm, [], stim, srate);
%figure;plot(t,stim,t,pred)

% Find the full width at half max value
halfMax = (min(pred) + max(pred)) / 2;
index1 = find(pred >= halfMax, 1, 'first');
index2 = find(pred >= halfMax, 1, 'last');
fwhmx = t(index2) - t(index1);

derivedPrm(3)    = fwhmx;

pred_transient = pred;


%%% OLD %%%
% paramNames = {'t2pk', 'rAsymp', 'rDouble', 'Tisi'};
% %% Compute derived parameter Rdouble
% 
% T = 2;
% stimLength = round(T*srate);
% stim_on = [1 2 4 8 16 32]/60;
% stim = zeros(stimLength,length(stim_on));
% 
% for ii = 1:length(stim_on)
%     stim_end = round(stim_on(ii)*srate);
%     stim(1:stim_end,ii) = 1;
% end    
%     
% [~, pred] = objFunction(prm, [], stim, srate);
% rsp       = sum(pred, 1);
% linpred   = rsp*2;
% rdub      = rsp(2:end)./linpred(1:end-1);
% 
% derivedPrm{3} = rdub;
% % % debug
% %figure;plot(stim, 'LineWidth', 2)
% %hold on; plot(pred, 'LineWidth', 2)
% 
% %% Compute derived parameter Tisi
% 
% T = 3;
% t_on = 100;
% stimLength = round(T*srate);
% stim = zeros(stimLength,1);
% stim(1 : t_on) = 1;
% 
% thresh = [0.8 0.9 0.95 1];
% 
% % compute summed response to first pulse
% [~, pred] = objFunction(prm, [], stim, srate);
% rsp = sum(pred(1:1000),1);
% 
% % create set of new stimuli with a second pulse 
% stim2 = repmat(stim,1,stimLength-2*t_on);
% pulse2_st = t_on+1:stimLength-t_on;
% for ii = 1:length(pulse2_st)
%     stim2(pulse2_st(ii) : pulse2_st(ii)+t_on-1, ii) = 1;
% end
% 
% % predict responses for 2 pulse stimuli
% [~, pred2] = objFunction(prm, [], stim2, srate);
% 
% % % debug
% % figure;
% % subplot(2,1,1);plot(pred2); yl = get(gca, 'YLim');
% % subplot(2,1,2);plot(pred2-pred); set(gca, 'YLim', yl);
% %  
% % subtract the response to first pulse
% pred2 = pred2-pred;
% % compute sum over matched windows relative to stimulus onset
% rsp2 = [];
% for ii = 1:size(pred2,2)-1000
%     rsp2(ii) = sum(pred2(pulse2_st(ii):pulse2_st(ii)+1000,ii),1);
% end
% %rsp2  = sum(pred2,1);
% 
% % compute t_isi
% isi_samples = nan(length(thresh),1);
% for ii = 1:length(thresh)
%     % find the 2 pulse condition for which the sum of the 2 pulse stimulus
%     % is equal to thresh * the first pulse
%     resp_diff = rsp2-(rsp*thresh(ii));
%     stim_inx = find(round(resp_diff)>=0);
%     if any(stim_inx) 
%     % compute t_isi (gap between offset of pulse 1 and onset of pulse 2)
%         isi_samples(ii) = (pulse2_st(stim_inx(1)) - t_on);
%     else
%         isi_samples(ii) = nan;
%     end
% end
% 
% t_isi = isi_samples./srate;
% derivedPrm{4} = t_isi;
% 
% % % debug
% stim_inx = isi_samples+t_on;
% resp_diff = rsp2-(rsp*thresh(1));

% figure('Position', [ 354    20   803   540]);%clf
% subplot(2,2,1);hold on
% plot(stim, 'LineWidth', 2)
% plot(pred, 'LineWidth', 2)
% set(gca, 'Xlim', [0 2000]);
% title('onepulse 100 ms');
% 
% subplot(2,2,2);hold on
% plot(resp_diff, 'LineWidth',2, 'Color','k');
% plot(stim_inx(1),resp_diff(stim_inx(1)), 'bo', 'LineWidth', 2, 'MarkerSize', 10);
% plot(stim_inx(2),resp_diff(stim_inx(2)), 'co', 'LineWidth', 2, 'MarkerSize', 10);
% plot(stim_inx(3),resp_diff(stim_inx(3)), 'mo', 'LineWidth', 2, 'MarkerSize', 10);
% plot(stim_inx(4),resp_diff(stim_inx(4)), 'ro', 'LineWidth', 2, 'MarkerSize', 10);
% legend({'difference', '80%', '90%', '95%', '100%'}, 'Location', 'SouthEast');
% title('sum(pulse 2) - sum(pulse 1)');
% set(gca, 'Xlim', [0 2000]);
% 
% stim_on  = find(stim == 1);
% stim_lth = length(stim_on);
% stim_end = stim_on(end);
% subplot(2,2,3);hold on
% pulse2_st  = stim_end + isi_samples(1) + 1;
% pulse2_end = stim_end + isi_samples(1) + 1 + stim_lth;
% stim2 = stim;
% stim2(pulse2_st : pulse2_end) = 1;
% [~, pred2] = objFunction(prm, [], stim2, srate);
% plot(stim2, 'LineWidth', 2)
% plot(pred2, 'LineWidth', 2)
% set(gca, 'Xlim', [0 2000]);
% title('80% recovery');
% 
% subplot(2,2,4);hold on
% pulse2_st  = stim_end + isi_samples(4) + 1;
% pulse2_end = stim_end + isi_samples(4) + 1 + stim_lth;
% stim2 = stim;
% stim2(pulse2_st : pulse2_end) = 1;
% [~, pred2] = objFunction(prm, [], stim2, srate);
% plot(stim2, 'LineWidth', 2)
% plot(pred2, 'LineWidth', 2)
% set(gca, 'Xlim', [0 2000]);
% title('100% recovery');

end


