
[nSamp, nStim, nDatasets] = size(data);

chan_inx = contains(channels.name, 'V1');

%% plot example of sub-additive temporal summation
figure('Position', [360    44   879   300]); hold on;

conditionsOfInterest = {'ONEPULSE-2', 'ONEPULSE-3'};
timepointsOfInterest = [-0.05 0.35];

stim_inx = contains(stim_info.name, conditionsOfInterest);
t_ind    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

s = stim_ts(t_ind,stim_inx);
d = data(t_ind,stim_inx,chan_inx);
sft = length(find(s(:,1)));
d_shift = padarray(d(:,1), [sft, 0], 0, 'pre');
s_shift = padarray(s(:,1), [sft, 0], 0, 'pre');
d_shift = d_shift(1:size(s,1));
s_shift = s_shift(1:size(s,1));
d_sum = sum([d_shift d(:,1)],2);
d_copy = d;
d_copy(:,1) = nan;
d_copy(:,2) = d_sum;
maxresp = max(d(:,1)); % scale stimulus to max of lowest duration

plot(flatten(s), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
plot(flatten(d./maxresp), 'k', 'LineWidth', 2);
plot(flatten(d_copy./maxresp), 'k--', 'LineWidth', 2);
set(gca, 'XTick',1:size(d,1):length(find(stim_inx))*size(d,1), 'XTickLabel', []);
set(gca, 'XTickLabel', stim_info.duration(stim_inx)); xlabel('duration (s)');
axis tight
title('sub-additive temporal summation');

%% plot example of adaptation
figure('Position', [360    44   879   300]); hold on;

conditionsOfInterest = {'TWOPULSE-2', 'TWOPULSE-6'};
timepointsOfInterest = [-0.05 1];

stim_inx = contains(stim_info.name, conditionsOfInterest);
t_ind    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

s = stim_ts(t_ind,stim_inx);
d = data(t_ind,stim_inx,chan_inx);
maxresp = max(d(:,1)); % scale stimulus to max of lowest duration

plot(flatten(s), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
plot(flatten(d./maxresp), 'k', 'LineWidth', 2);
set(gca, 'XTick',1:size(d,1):length(find(stim_inx))*size(d,1), 'XTickLabel', []);
set(gca, 'XTickLabel', stim_info.ISI(stim_inx)); xlabel('ISI (s)');
axis tight
title('adaptation to repeated stimulus');

%% plot subset of contrast 
figure('Position', [360    44   879   300]); hold on;

stim_inx = contains(stim_info.name, {'CRF-1','CRF-5'});
s = stim_ts(:,stim_inx);
d = data(:,stim_inx,chan_inx);
maxresp = max(d(:)); % scale stimulus to max of all included conditions

plot(flatten(s*maxresp), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
plot(flatten(d), 'k', 'LineWidth', 2);
set(gca, 'XTick',1:size(d,1):length(find(stim_inx))*size(d,1), 'XTickLabel', []);
set(gca, 'XTickLabel', stim_info.contrast(stim_inx)*100); xlabel('contrast (%)');
axis tight
title('different dynamics with low contrast');

%% plot contrast conditions on top of each other and normalized
figure('Position', [360    44   879   300]); hold on;

conditionsOfInterest = {'CRF-1','CRF-2', 'CRF-3', 'CRF-4','CRF-5'};
%conditionsOfInterest = {'CRF-1','CRF-5'};

timepointsOfInterest = [-0.05 0.25];


stim_inx = contains(stim_info.name, conditionsOfInterest);
t_ind    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

s = stim_ts(t_ind,stim_inx);
d = data(t_ind,stim_inx,chan_inx);
maxresp = max(d(:)); % scale stimulus to max of lowest duration
subplot(1,2,1);
colors = flipud(gray(length(conditionsOfInterest)+1)); colors = colors(2:end,:);
ecog_plotMultipleTimeCourses(t(t_ind), d, [], colors);
ecog_plotMultipleTimeCourses(t(t_ind), s*maxresp, [], colors, [], [], [-2 maxresp]);
legend(num2str(stim_info.contrast(stim_inx)*100));
title('contrast amplitude shift');

d = d./max(d);
subplot(1,2,2);
colors = flipud(gray(length(conditionsOfInterest)+1)); colors = colors(2:end,:);
ecog_plotMultipleTimeCourses(t(t_ind), d, [], colors, [], [], [-0.1 1]);
plot(t(t_ind), s(:,end), 'k','LineWidth', 2)
legend(num2str(stim_info.contrast(stim_inx)));
title('contrast latency shift');

%% contrast conditions normalized with stimulus timecourse
figure; hold on;

stim_inx = contains(stim_info.name, {'CRF-5'});
s = stim_ts(:,stim_inx);
d = data(:,stim_inx,chan_inx); d = d./max(d);
plot(t,s, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
plot(t,d, 'k', 'LineWidth', 2);

stim_inx = contains(stim_info.name, {'CRF-1'});
d = data(:,stim_inx,chan_inx); d = d./max(d);
maxresp = max(d(:)); % scale stimulus to max of all included conditions
plot(t,d, 'Color', [0.8 0.8 0.8], 'LineWidth', 2);
