% Figure 4 : 
% Contrast-dependent temporal dynamics of neuronal responses in human V1.
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
stim = D.stim;
stim_info = D.stim_info;

% Prepare figure
figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

% Subplot positions: % [left bottom width height]
posb = [0.05 0.1 0.92 0.3];

posa1 = [0.04 0.5 0.15 0.45];
posa2 = [0.22 0.5 0.15 0.45];
posc1 = [0.43 0.5 0.15 0.45];
posc2 = [0.63 0.5 0.15 0.45];
posc3 = [0.82 0.5 0.15 0.45];

%% Panel A: contrast response timecourses (original and normalized)

% Select conditions to plot
conditionsOfInterest = {'CRF-1','CRF-2', 'CRF-3', 'CRF-4','CRF-5'}; 

% Select time window to plot
timepointsOfInterest = [-0.05 0.31];

% Look up corresponding indices in the data
stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

% Generate time courses + linear prediction
d = data(t_idx,stim_idx,1);

% Make legend
x = stim_info.contrast(stim_idx)*100; %  in percent
l = [];
for ii = length(x):-1:1
    if ii == 1,  l{length(x)+1-ii} = sprintf('%-5.2f', x(ii)); end
    if ii == 2,  l{length(x)+1-ii} = sprintf('%-5.1f', x(ii)); end
    if ii > 2,   l{length(x)+1-ii} = sprintf('%-5.0f', x(ii)); end
end

% Plot
subplot('position', posa1); cla; hold on
colors = gray(length(conditionsOfInterest)+1); colors = colors(1:end-1,:);
ecog_plotMultipleTimeCourses(t(t_idx)*1000, fliplr(d), [], colors);
set(gca, 'xtick', [0 300]);
set(gca, 'ylim', [-2 22], 'ytick', []);
xlabel('Time (ms)'); ylabel('Neural response');

legend(l);
legend('boxoff')

d = d./max(d);
subplot('position', posa2); hold on
ecog_plotMultipleTimeCourses(t(t_idx)*1000, fliplr(d), [], colors, [], [], [-0.1 1]);
set(gca, 'xtick', [0 300]);
set(gca, 'ylim', [-0.1 1.15]);
set(gca, 'ytick', [0 1], 'ytick', []);
xlabel('Time (ms)'); ylabel('Neural response (normalized)');

%% Panel B: data and fits

conditionsOfInterest = {'CRF'};
timepointsOfInterest = [-0.1 1];

stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

s = stim(t_idx,stim_idx);
d = data(t_idx,stim_idx);
d_se_l = data_se(t_idx,stim_idx,1);
d_se_u = data_se(t_idx,stim_idx,2);
p = pred(t_idx,stim_idx);
p_se_l = pred_se(t_idx,stim_idx,1);
p_se_u = pred_se(t_idx,stim_idx,2);

maxresp = max(d(:)); % scale stimulus to max 

% Plot
subplot('position', posb); cla; hold on
hs = plot(s(:)*maxresp, 'color', [0.7 0.7 0.7], 'linewidth', 1);
hcid = ciplot(d_se_l(:), d_se_u(:), [], 'k', 0.25);
hd = plot(d(:), 'k-', 'linewidth', 2);
hcip = ciplot(p_se_l(:), p_se_u(:), [], 'r', 0.25);
hp = plot(p(:), 'r-', 'linewidth', 2);
set(get(get(hcid,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(hcip,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xtick',1:size(d,1):size(d,2)*size(d,1), 'xticklabel', stim_info.contrast(stim_idx)*100);
box off,  axis tight

set(gca, 'ylim', [-2 20]);
set(gca, 'xlim', [-20 length(s(:)) + 20]);

xlabel('Contrast (%)'); ylabel('Change in power (x-fold)'); 
legend({'Stimulus','Neural response', 'DN prediction'}, 'location', 'northwest');
legend('boxoff');

%% Panel C: contrast dynamics quantified across electrodes

% Find stimulus index
stim_idx = find(contains(stim_info.name, 'CRF'));
x = stim_info.contrast(stim_idx) * 100; 

% Generate new stimuli based on params
nStim = 50;
[stim2, stim_info2] = tde_simulateNewStimuli(t,nStim);
stim_idx2 = find(contains(stim_info2.name, 'CRF'));
x2 = stim_info2.contrast(stim_idx2) * 100; 

% Predict model responses for new stimuli
srate = D.channels.sampling_frequency(1);
pred = nan([size(stim2(:,stim_idx2)) height(D.channels)]);
for ii = 1:size(D.params,2)
    prm = D.params(:,ii);
   	[~, pred(:,:,ii)] = D.objFunction(prm, [], stim2(:,stim_idx2), srate);      
end

% Determine time index over which to compute summary statistics
t_idx = t>0.05 & t<1.0;

% Concatenate data and prediction (in order to make sure probabilistic
% assignment of electrodes to areas is done the same way for both).
d = D.data(t_idx,stim_idx,:);
d = cat(2,d,pred(t_idx, :,:));

% Compute sum across stim_on window
m_conc = squeeze(sum(d,1)); 
[m_conc, se_conc] = averageWithinArea(m_conc, group_prob, @median, numboot);

% Plot

% 1. Contrast response function
subplot('position', posc1); cla; hold on

% Plot linear prediction 
h0 = line([0 x(end)], [0 1], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Plot data
nStim = length(stim_idx);
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'errbar', 1, [], 50);

% Plot prediction
m = m_conc(nStim+1:end);
se = se_conc(nStim+1:end,:);
tde_plotPoints(m, se, x2, 'ci', 1, [], 50,'r');

% Format axes
ticklabelsX = num2str(x); ticklabelsX(2:4,:) = ' ';
set(gca, 'xlim', [0 105], 'xtick', x, 'xticklabel', ticklabelsX);
xlabel('Contrast (%)'); ylabel('Summed broadband power (0-1s)'); 

% 2. Time to peak
subplot('position', posc2); cla; hold on

% Compute peak across trial window
[~,I] = max(d,[],1);
T = t(t_idx);
[m_conc, se_conc] = averageWithinArea(squeeze(T(I)), group_prob,  @median, numboot);

m_conc = m_conc * 1000; % convert to ms
se_conc = se_conc * 1000;

% Plot data
nStim = length(stim_idx);
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'errbar', 0, [], 50);

% Plot prediction
m = m_conc(nStim+1:end);
se = se_conc(nStim+1:end,:);
tde_plotPoints(m, se, x2, 'ci', 0, [], 50, 'r');

% Format axes
l = get(gca, 'YLim'); ylim([50 l(2)]);
ticklabelsX = num2str(x); ticklabelsX(2:4,:) = ' ';
set(gca, 'xlim', [0 105], 'xtick', x, 'xticklabel', ticklabelsX);
xlabel('Contrast (%)'); ylabel('Time-to-peak (ms)'); 

% 3. Ratio of sustained to transient
subplot('position', posc3); cla; hold on
T = D.t(t_idx);

% Smooth time courses to get better estimates of max and offset response levels
% Apply same amount of smoothing to both data and model predictions
for ii = 1:size(d,2)
    for jj = 1:size(d,3)
        d(:,ii,jj) = smooth(d(:,ii,jj),150);
    end
end

% DEBUG
%figure; plot(T,squeeze(d(:,5,group_prob>0)));

[M] = max(d,[],1); % value at peak
t_off = (T == 0.5); % offset timepoint
O = d(t_off,:,:); % value at offset
R = squeeze(O./M); % divide value at offset with value at peak

[m_conc, se_conc] = averageWithinArea(R, group_prob, @median, numboot);

% Plot data
nStim = length(stim_idx);
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'errbar', 0, [], 50);

% Plot prediction
m = m_conc(nStim+1:end);
se = se_conc(nStim+1:end,:);
tde_plotPoints(m, se, x2, 'ci', 0, [], 50, 'r');

% Format axes
ticklabelsX = num2str(x); ticklabelsX(2:4,:) = ' ';
set(gca, 'xlim', [0 105], 'xtick', x, 'xticklabel', ticklabelsX);
set(gca, 'ylim', [0 1]);
xlabel('Contrast (%)'); ylabel('Ratio offset/peak'); 

set(findall(gcf,'-property','FontSize'),'FontSize',24)
