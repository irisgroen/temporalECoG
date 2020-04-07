function [derivedPrm, pred, t] = tde_computeDerivedParams(objFunction, prm)
% Simulate model response to a long sustained stimulus
% Extract two summary statistics:
% - t2pk = time to peak 
% - r_asymp = magnitude or response at asymptote

%% Compute derived parameters Time2Peak and Rasymptote

t    = 0.001 : 0.001 : 10;
stim = ones(length(t),1);
srate = 1/median(diff(t)); % samples per second

[~, pred] = objFunction(prm, [], stim, srate);

[~,x] = max(pred);
derivedPrm.t2pk    = t(x);
derivedPrm.r_asymp = pred(end)/max(pred);

%% Compute derived parameter Rdouble

T = 1;
stimLength = round(T*srate);
stim = zeros(stimLength,2);
stim(1:100,1) = 1;
stim(1:200,2) = 1;

[~, pred] = objFunction(prm, [], stim, srate);
rsp = sum(pred, 1);
derivedPrm.r_double = rsp(2)/(2*rsp(1));

% % debug
%figure;plot(stim, 'LineWidth', 2)
%hold on; plot(pred, 'LineWidth', 2)

%% Compute derived parameter Tisi

T = 3;
t_on = 100;
stimLength = round(T*srate);
stim = zeros(stimLength,1);
stim(1 : t_on) = 1;

thresh_full = 1; % 90 percent recovered
thresh_part = 0.90; % 90 percent recovered

% compute summed response to first pulse
[~, pred] = objFunction(prm, [], stim, srate);
rsp = sum(pred,1);

% create set of new stimuli with a second pulse 
stim2 = repmat(stim,1,stimLength-2*t_on);
pulse2_st = t_on+1:stimLength-t_on;
for ii = 1:length(pulse2_st)
    stim2(pulse2_st(ii) : pulse2_st(ii)+t_on-1, ii) = 1;
end

% predict responses for 2 pulse stimuli
[~, pred2] = objFunction(prm, [], stim2, srate);
rsp2 = sum(pred2,1);

% find the 2 pulse condition for which the sum of the 2 pulse stimulus is
% equal to first pulse plus thresh * the first pulse
resp_diff_full = rsp2-(rsp+rsp*thresh_full);
resp_diff_part = rsp2-(rsp+rsp*thresh_part);
stim_inx_full = find(round(resp_diff_full)>=0);
stim_inx_part = find(round(resp_diff_part)>=0);

% compute t_isi (gap between offset of pulse 1 and onset of pulse 2)
isi_samples_full = (pulse2_st(stim_inx_full(1)) - t_on);
isi_samples_part = (pulse2_st(stim_inx_part(1)) - t_on);
t_isi = isi_samples_full/srate;
derivedPrm.t_isi = t_isi;

% % debug
figure('Position', [ 354    20   803   540]);%clf
subplot(2,2,1);hold on
plot(stim, 'LineWidth', 2)
plot(pred, 'LineWidth', 2)
set(gca, 'Xlim', [0 2000]);
title('onepulse 100 ms');

subplot(2,2,2);hold on
plot(resp_diff_full, 'LineWidth',2, 'Color','k');
plot(stim_inx_full(1),resp_diff_full(stim_inx_full(1)), 'bo', 'LineWidth', 2, 'MarkerSize', 10);
plot(stim_inx_part(1),resp_diff_full(stim_inx_part(1)), 'co', 'LineWidth', 2, 'MarkerSize', 10);
legend({'difference', '100%', '90%'}, 'Location', 'SouthEast');
title('twopulse-(onepulse+thresh*onepulse)');
set(gca, 'Xlim', [0 2000]);

stim_on  = find(stim == 1);
stim_lth = length(stim_on);
stim_end = stim_on(end);
subplot(2,2,3);hold on
pulse2_st  = stim_end + isi_samples_part + 1;
pulse2_end = stim_end + isi_samples_part + 1 + stim_lth;
stim2 = stim;
stim2(pulse2_st : pulse2_end) = 1;
[~, pred2] = objFunction(prm, [], stim2, srate);
plot(stim2, 'LineWidth', 2)
plot(pred2, 'LineWidth', 2)
set(gca, 'Xlim', [0 2000]);
title('90% recovery');

subplot(2,2,4);hold on
pulse2_st  = stim_end + isi_samples_full + 1;
pulse2_end = stim_end + isi_samples_full + 1 + stim_lth;
stim2 = stim;
stim2(pulse2_st : pulse2_end) = 1;
[~, pred2] = objFunction(prm, [], stim2, srate);
plot(stim2, 'LineWidth', 2)
plot(pred2, 'LineWidth', 2)
set(gca, 'Xlim', [0 2000]);
title('100% recovery');

end


