% tde_mkFigure2

% Load data and fits

% electrode-averaged data and DN model fits
modelfun = @LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datatype = 'electrodeaverages';
[d1] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% individual electrodes and DN model fits
datatype = 'individualelecs';
[d2] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% Load stimulus info
[stim_ts, stim_info] = tde_generateStimulusTimecourses(d1.options.stimnames,d1.t);
stim_info.duration = stim_info.duration*1000; % convert to ms

% Subplot positions: % [left bottom width height]
posa = [0.1 0.55 0.4 0.3];
posb = [0.1 0.1 0.8 0.25];
posc = [0.55 0.45 0.35 0.45];

figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

%% Panel A: example of compressive temporal summation 

% Select two conditions to plot
conditionsOfInterest = {'ONEPULSE-1', 'ONEPULSE-2'};

% Select time window to plot
timepointsOfInterest = [-0.05 0.3];

% Look up corresponding indices in the data
stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = d1.t>timepointsOfInterest(1) & d1.t<=timepointsOfInterest(2);

% Generate time courses + linear prediction
s = stim_ts(t_idx,stim_idx);
d = d1.data(t_idx,stim_idx,1);
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

% For plotting purposes, add a little space after the first stimulus
dummy = nan(40,2);
s = cat(1, s, dummy);
d = cat(1, d, dummy);
d_copy = cat(1, d_copy, dummy);

% Plot
subplot('position', posa); hold on
plot(s(:), 'color', [0.5 0.5 0.5], 'lineWidth', 1);
plot((d(:)./maxresp), 'k', 'lineWidth', 2);
plot((d_copy(:)./maxresp), 'k:', 'lineWidth', 2);

% Set axes
set(gca, 'xtick',1:size(d,1):length(find(stim_idx))*size(d,1), 'ytick', []);
set(gca, 'xticklabel', stim_info.duration(stim_idx)); box off
xlabel('Stimulus duration (ms)'); ylabel('Neural response'); title('Temporal summation is sub-additive', 'fontsize', 20); 

% Add legend
legend({'Stimulus', 'Neural data', 'Linear prediction'}, 'location', 'northwest', 'fontsize', 18);
legend('boxoff')

%% Panel B: data and fits
conditionsOfInterest = {'ONEPULSE'};
timepointsOfInterest = [-0.1 0.8];

stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = d1.t>timepointsOfInterest(1) & d1.t<=timepointsOfInterest(2);

s = stim_ts(t_idx,stim_idx);
d = d1.data(t_idx,stim_idx,1);
p = d1.pred(t_idx,stim_idx,1);
maxresp = max(d(:,1)); % scale stimulus to max of first condition

subplot('position', posb); hold on
hs = plot(s(:)*maxresp, 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(d(:), 'k-', 'linewidth', 2);
hp = plot(p(:), 'r-', 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xtick',1:size(d,1):size(d,2)*size(d,1), 'xticklabel', stim_info.duration(stim_idx));
box off,  axis tight

ylim([-2 20]);
xlim([-20 length(s(:)) + 20]);

xlabel('Stimulus duration (ms)'); ylabel('Change in power (x-fold)'); title('Broadband responses to increasing durations', 'fontsize', 20); 
legend({'Neural data', 'DN model prediction'}, 'location', 'northwest', 'fontsize', 18);
legend('boxoff');

%% Panel C: temporal summation across electrodes

% Find stimulus index
stim_idx = find(contains(stim_info.name, 'ONEPULSE'));
x = stim_info.duration(stim_idx); % in ms

% Generate new stimuli based on params
nStim = 50;
[stim2, stim_info2] = tde_simulateNewStimuli(d2.t,nStim);
stim_info2.duration = stim_info2.duration*1000; % convert to ms
stim_idx2 = find(contains(stim_info2.name, 'ONEPULSE'));
x2 = stim_info2.duration(stim_idx2); % in ms

% Predict model responses for new stimuli
srate = d2.channels.sampling_frequency(1);
pred = nan([size(stim2(:,stim_idx2)) height(d2.channels)]);
for ii = 1:size(d2.params,2)
    prm = d2.params(:,ii);
   	[~, pred(:,:,ii)] = d2.objFunction(prm, [], stim2(:,stim_idx2), srate);      
end

% Determine time index over which to compute summary statistics
t_idx = d2.t>0 & d2.t<1;

% Concatenate data and prediction (in order to make sure probabilistic
% assignment of electrodes to areas is done the same way for both).
d = d2.data(t_idx,stim_idx,:);
d = cat(2,d,pred(t_idx, :,:));

% Compute sum across stim_on window
m_conc = squeeze(sum(d,1)); 
[~, channels, group_prob] = groupElecsByVisualArea(d2.channels, 'probabilisticresample', {'V1'});   
[m_conc, se_conc] = averageWithinArea(m_conc, group_prob, [], 10000);

% Plot
subplot('position', posc); hold on

% Plot linear prediction 
h0 = line([0 x(end)], [0 1], 'linestyle', ':', 'LineWidth', 2, 'color', [0 0 0]);
set(gca, 'xtick', x);

% Plot data
nStim = length(stim_idx);
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'errbar', 1)

% Plot prediction
m = m_conc(nStim+1:end);
se = se_conc(nStim+1:end,:);
tde_plotPoints(m, se, x2, 'ci', 1)

% Format axes
xlabel('Stimulus duration (ms)'); ylabel('Summed broadband timecourse (0-1s)'); title('Temporal summation', 'fontsize', 20); 
legend({'Linear prediction', 'Neural data', 'DN model prediction'}, 'location', 'southeast', 'fontsize', 18);

legend('boxoff')
axis square
axis tight

set(findall(gcf,'-property','FontSize'),'FontSize',20)

