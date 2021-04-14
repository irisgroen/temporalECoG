function [derivedPrm, paramNames, pred_sustained, pred_transient, t] = tde_computeDerivedParams(objFunction, prm)
% Simulate model response to a long sustained stimulus
% Extract four summary statistics:
% - t2pk    = time to peak 
% - r_asymp = magnitude or response at asymptote
% - Rdouble = changed in summed response with doubling of duration
% - t_isi = magnitude or response at asymptote

paramNames = {'t2pk', 'rAsymp', 'wSize', 'Tisi', 'c50'};
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

%% Compute derived parameters Tisi (time to achieve full recovery)

T = 2;
t_on = 134;
stimLength = round(T*srate);
stim = zeros(stimLength,1);
stim(1 : t_on) = 1;

thresh = [1];

% compute summed response to first pulse
[~, pred] = objFunction(prm, [], stim, srate);
rsp = max(pred(1:1000));

% create set of new stimuli with a second pulse 

stim2 = repmat(stim,1,(stimLength-2*t_on)/2);
pulse2_st = t_on+1:2:stimLength-t_on;
for ii = 1:length(pulse2_st)
    stim2(pulse2_st(ii) : pulse2_st(ii)+t_on-1, ii) = 1;
end

% predict responses for 2 pulse stimuli
[~, pred2] = objFunction(prm, [], stim2, srate);

% % debug
% figure;
% subplot(2,1,1);plot(pred2); yl = get(gca, 'YLim');
% subplot(2,1,2);plot(pred2-pred); set(gca, 'YLim', yl);
%  
% subtract the response to first pulse
pred2 = pred2-pred;

% compute max over matched windows relative to stimulus onset
rsp2 = [];
for ii = 1:size(pred2,2)-500
    rsp2(ii) = max(pred2(pulse2_st(ii):pulse2_st(ii)+1000,ii));
end

% compute t_isi
isi_samples = nan(length(thresh),1);
for ii = 1:length(thresh)
    % find the 2 pulse condition for which the sum of the 2 pulse stimulus
    % is equal to thresh * the first pulse
    resp_diff = rsp2-(rsp*thresh(ii));
    stim_inx = find(round(resp_diff)>=0);
    if any(stim_inx) 
    % compute t_isi (gap between offset of pulse 1 and onset of pulse 2)
        isi_samples(ii) = (pulse2_st(stim_inx(1)) - t_on);
    else
        isi_samples(ii) = nan;
    end
end

t_isi = isi_samples./srate;
derivedPrm(4) = t_isi;

%% Compute derived parameter C50 (contrast to reach 50% output)

% Simulate response to various contrast levels
t    = 0.001 : 0.001 : 0.5;
srate = 1/median(diff(t)); % samples per second

% Generate stimuli with different contrast levels
stim = nan(length(t),100);
for ii = 1:100
    stim(:,ii) = ones(length(t),1)*ii*0.01;
end

% Predict responses
[~, pred] = objFunction(prm, [], stim, srate);

% Find maximum response at each contrast level
maxResp = max(pred);
% Find 50% of the maximum at the 100% contrast level
halfMax = 0.5*maxResp(:,100);
% Find the contrast level corresponding to 50% of max
resp_diff = maxResp-halfMax;
derivedPrm(5) = find(resp_diff> 0,1);

%% Compute derived parameter Cscale (to what extent predictions match just by scaling)

end


