%tde_mkFigure 4

% Load data and fits

% electrode-averaged data and DN model fits
modelfun = @LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datatype = 'electrodeaverages';
[data1] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% individual electrodes and DN model fits
datatype = 'individualelecs';
[data2] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

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
