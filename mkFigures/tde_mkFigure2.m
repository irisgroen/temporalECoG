% Figure 2 : 
% Sub-additive temporal summation in neural responses in human V1.
%
% 2022 Iris Groen

% Load data and fits

modelfun = @DN;
xvalmode = 0;
datatype = 'individualelecs';
[D] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% Select electrodes and compute averages 
numboot = 10000;
[~, ~, group_prob] = groupElecsByVisualArea(D.channels, 'probabilisticresample', {'V1'});   
[data, data_se] = averageWithinArea(D.data, group_prob, @mean, numboot);
[pred, pred_se] = averageWithinArea(D.pred, group_prob, @mean, numboot);
t = D.t;
stim_ts = D.stim;
stim_info = D.stim_info;

% Prepare figure
figure(4); clf
set(gcf,'position',get(0, 'screensize'));

% Subplot positions: % [left bottom width height]
posa = [0.1 0.55 0.35 0.3];
posb = [0.55 0.5 0.35 0.45];
posc = [0.1 0.1 0.8 0.25];

%% Panel A: example of compressive temporal summation 

% Select two conditions to plot
conditionsOfInterest = {'ONEPULSE-1', 'ONEPULSE-2'};

% Select time window to plot
timepointsOfInterest = [-0.05 0.3];

% Look up corresponding indices in the data
stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

% Generate time courses + linear prediction
s = stim_ts(t_idx,stim_idx);
d = data(t_idx,stim_idx,1);
sft = length(find(s(:,1)));
d_shift = padarray(d(:,1), [sft, 0], 0, 'pre');
s_shift = padarray(s(:,1), [sft, 0], 0, 'pre');
d_shift = d_shift(1:size(s,1));
s_shift = s_shift(1:size(s,1));
d_sum = sum([d_shift d(:,1)],2);
d_copy = d;
d_copy(:,1) = nan;
d_copy(:,2) = d_sum;
maxresp = max(d(:,1)); % Scale stimulus to max of lowest duration

% For plotting purposes, add a little space after the first stimulus
dummy = nan(40,2);
s = cat(1, s, dummy);
d = cat(1, d, dummy);
d_copy = cat(1, d_copy, dummy);

% Plot
subplot('position', posa); cla; hold on
plot(s(:), 'color', [0.5 0.5 0.5], 'lineWidth', 1);
plot((d(:)./maxresp), 'k', 'lineWidth', 2);
plot((d_copy(:)./maxresp), 'k:', 'lineWidth', 2);

% Set axes
set(gca, 'xtick',1:size(d,1):length(find(stim_idx))*size(d,1), 'ytick', []);
set(gca, 'xticklabel', stim_info.duration(stim_idx)); box off
xlabel('Stimulus duration (ms)'); ylabel('Response magnitude');

% Add legend
legend({'Stimulus', 'Neural response', 'Linear prediction'}, 'location', 'northwest');
legend('boxoff')

%% Panel B: temporal summation across electrodes

% Find stimulus index
stim_idx = find(contains(stim_info.name, 'ONEPULSE'));
x = stim_info.duration(stim_idx)*1000; % in ms

% Generate new stimuli based on params
nStim = 50;
[stim2, stim_info2] = tde_simulateNewStimuli(t,nStim);
stim_info2.duration = stim_info2.duration*1000; % convert to ms
stim_idx2 = find(contains(stim_info2.name, 'ONEPULSE'));
x2 = stim_info2.duration(stim_idx2); % in ms

% Predict model responses for new stimuli
srate = D.channels.sampling_frequency(1);
pred2 = nan([size(stim2(:,stim_idx2)) height(D.channels)]);
for ii = 1:size(D.params,2)
    prm = D.params(:,ii);
   	[~, pred2(:,:,ii)] = D.objFunction(prm, [], stim2(:,stim_idx2), srate);      
end

% Determine time index over which to compute summary statistics
t_idx = t>0 & t<1;

% Concatenate data and prediction (in order to make sure probabilistic
% assignment of electrodes to areas is done the same way for both).
d = D.data(t_idx,stim_idx,:);
d = cat(2,d,pred2(t_idx, :,:));

% Compute sum across stim_on window
m_conc = squeeze(sum(d,1)); 
[m_conc, se_conc] = averageWithinArea(m_conc, group_prob, @median, numboot);

% Plot
subplot('position', posb);cla; hold on

% Calculate data
nStim = length(stim_idx);
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);

% Calculate linear prediction 
%lin_pred = [m(1) m(1)*2 m(1)*4 m(1)*8 m(1)*16 m(1)*32];
lin_pred = [0 1];

% Plot linear prediction 
h0 = line([0 x(end)], [0 lin_pred(end)], 'linestyle', ':', 'LineWidth', 2, 'color', [0 0 0]);

% Plot data
tde_plotPoints(m, se, x, 'errbar', 1,[],50);

% Plot prediction
m = m_conc(nStim+1:end);
se = se_conc(nStim+1:end,:);
tde_plotPoints(m, se, x2, 'ci', 1,[],50, 'r');

% Format axes
xlabel('Stimulus duration (ms)'); ylabel('Summed broadband power (0-1s)'); 
legend({'Linear prediction', 'Neural response', 'DN prediction'}, 'location', 'southeast');
set(gca, 'xtick', x, 'xticklabelrotation', 45);

legend('boxoff')
axis tight
xlim([0 550]);
ylim([0 1.3]);
axis square

%% Panel C: data and fits

% Select stimuli and time-points of interest
conditionsOfInterest = {'ONEPULSE'};
timepointsOfInterest = [-0.1 0.8];

stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

% Generate time courses and predictions
s = stim_ts(t_idx,stim_idx);
d = data(t_idx,stim_idx);
d_se_l = data_se(t_idx,stim_idx,1);
d_se_u = data_se(t_idx,stim_idx,2);
p = pred(t_idx,stim_idx);
p_se_l = pred_se(t_idx,stim_idx,1);
p_se_u = pred_se(t_idx,stim_idx,2);
maxresp = max(d(:,1)); % Scale stimulus to max of first condition

% Plot
subplot('position', posc); cla; hold on

hs = plot(s(:)*maxresp, 'color', [0.7 0.7 0.7], 'linewidth', 1);
hcid = ciplot(d_se_l(:), d_se_u(:), [], 'k', 0.25);
hd = plot(d(:), 'k-', 'linewidth', 2);
hcip = ciplot(p_se_l(:), p_se_u(:), [], 'r', 0.25);
hp = plot(p(:), 'r-', 'linewidth', 2);
set(get(get(hcid,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(hcip,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Format axes
set(gca, 'xtick',1:size(d,1):size(d,2)*size(d,1), 'xticklabel', stim_info.duration(stim_idx)*1000, 'xticklabelrotation', 45);
box off,  axis tight

ylim([-2 25]);
xlim([-20 length(s(:)) + 20]);

xlabel('Stimulus duration(ms)'); ylabel('Change in power (x-fold)');
legend({'Stimulus', 'Neural response', 'DN prediction'}, 'location', 'northwest');
legend('boxoff');

set(findall(gcf,'-property','FontSize'),'FontSize',24)

