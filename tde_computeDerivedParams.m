function [derivedPrm, paramNames, pred_sustained, pred_transient, t] = tde_computeDerivedParams(objFunction, prm)
% Simulate model response to a long sustained stimulus
% Extract four summary statistics:
% - t2pk    = time to peak 
% - r_asymp = magnitude or response at asymptote
% - Rdouble = changed in summed response with doubling of duration
% - t_isi = magnitude or response at asymptote

paramNames = {'t2pk', 'rAsymp', 'wSize', 'c50', 'cScale', 'cScaleR2'};
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

% Simulate response to a short (17 ms) pulse 
t    = 0.001 : 0.001 : 0.5;
srate = 1/median(diff(t)); % samples per second
stim = zeros(length(t),1);
stim(1:16) = 1;

[~, pred] = objFunction(prm, [], stim, srate);
%figure;plot(t,stim,t,pred)

% Find the full width at half max value
halfMax = (min(pred) + max(pred)) / 2;
index1 = find(pred >= halfMax, 1, 'first');
index2 = find(pred >= halfMax, 1, 'last');
fwhmx = t(index2) - t(index1);

derivedPrm(3)    = fwhmx;

pred_transient = pred;

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
derivedPrm(4) = find(resp_diff> 0,1);

%% Compute derived parameter Cscale (to what extent predictions match just by scaling)
X = pred(:,7); % 7% contrast
y = pred(:,100); % 100% contrast
stats = regstats(y,X);

% % debug
% figure;plot(t,X,t,y,t,stats.r)
% legend('low','high','resid')
% title(stats.mse);
derivedPrm(5) = stats.mse;
derivedPrm(6) = stats.rsquare;

end


