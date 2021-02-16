%tde_mkFigure 4

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

% Subplot positions: % [left bottom width height]
posb = [0.05 0.1 0.9 0.25];

posa1 = [0.04 0.45 0.15 0.45];
posa2 = [0.21 0.45 0.15 0.45];
posc1 = [0.45 0.45 0.15 0.45];
posc2 = [0.63 0.45 0.15 0.45];
posc3 = [0.80 0.45 0.15 0.45];

figure(1); clf
%figure;hold on
set(gcf, 'position',  get(0, 'screensize'));

%% Panel A: contrast response timecourses (original and normalized)

% Select conditions to plot
conditionsOfInterest = {'CRF-1','CRF-2', 'CRF-3', 'CRF-4','CRF-5'}; 
%conditionsOfInterest = {'CRF-2','CRF-3'};

% Select time window to plot
timepointsOfInterest = [-0.05 0.31];

% Look up corresponding indices in the data
stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = d1.t>timepointsOfInterest(1) & d1.t<=timepointsOfInterest(2);

% Generate time courses + linear prediction
s = d1.stim(t_idx,stim_idx);
d = d1.data(t_idx,stim_idx,1);
maxresp = max(d(:)); % scale stimulus to max of lowest duration

% Make legend
x = stim_info.contrast(stim_idx)*100; %  in percent
l = [];
for ii = length(x):-1:1
    if ii == 1,  l{length(x)+1-ii} = sprintf('%-5.2f', x(ii)); end
    if ii == 2,  l{length(x)+1-ii} = sprintf('%-5.1f', x(ii)); end
    if ii > 2,   l{length(x)+1-ii} = sprintf('%-5.0f', x(ii)); end
    %if ii == 1,  l{length(x)+1-ii} = sprintf('%-5.1f', x(ii)); end
    %if ii > 1,   l{length(x)+1-ii} = sprintf('%-5.0f', x(ii)); end
end

% Plot
subplot('position', posa1); hold on
colors = gray(length(conditionsOfInterest)+1); colors = colors(1:end-1,:);
ecog_plotMultipleTimeCourses(d1.t(t_idx)*1000, fliplr(d), [], colors);
set(gca, 'xtick', [0 300]);
set(gca, 'ylim', [-2 22], 'ytick', []);
xlabel('Time (ms)'); ylabel('Neural response');
title('Peak amplitude decrease','fontsize', 18);

legend(l,  'fontsize', 18);
legend('boxoff')


d = d./max(d);
subplot('position', posa2); hold on
ecog_plotMultipleTimeCourses(d1.t(t_idx)*1000, fliplr(d), [], colors, [], [], [-0.1 1]);
set(gca, 'xtick', [0 300]);
set(gca, 'ylim', [-0.1 1.15]);
set(gca, 'ytick', [0 1], 'ytick', []);
xlabel('Time (ms)'); ylabel('Neural response (normalized to peak)');
title('Peak latency shift', 'fontsize', 18);

%% Panel B: data and fits

conditionsOfInterest = {'CRF'};
timepointsOfInterest = [-0.1 1];

stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = d1.t>timepointsOfInterest(1) & d1.t<=timepointsOfInterest(2);

s = stim_ts(t_idx,stim_idx);
d = d1.data(t_idx,stim_idx,1);
p = d1.pred(t_idx,stim_idx,1);
maxresp = max(d(:)); % scale stimulus to max 

subplot('position', posb); hold on
hs = plot(s(:)*maxresp, 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(d(:), 'k-', 'linewidth', 2);
hp = plot(p(:), 'r-', 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xtick',1:size(d,1):size(d,2)*size(d,1), 'xticklabel', stim_info.contrast(stim_idx)*100);
box off,  axis tight

set(gca, 'ylim', [-2 20]);
set(gca, 'xlim', [-20 length(s(:)) + 20]);

xlabel('Stimulus contrast (%)'); ylabel('Change in power (x-fold)'); 
title('Broadband responses to increasing contrast', 'fontsize', 18); 
legend({'Neural data', 'DN model prediction'}, 'location', 'northwest', 'fontsize', 18);
legend('boxoff');

%% Panel C: contrast dynamics quantified across electrodes

% Find stimulus index
stim_idx = find(contains(stim_info.name, 'CRF'));
x = stim_info.contrast(stim_idx) * 100; 

% Generate new stimuli based on params
nStim = 50;
[stim2, stim_info2] = tde_simulateNewStimuli(d2.t,nStim);
stim_idx2 = find(contains(stim_info2.name, 'CRF'));
x2 = stim_info2.contrast(stim_idx2) * 100; 

% Predict model responses for new stimuli
srate = d2.channels.sampling_frequency(1);
pred = nan([size(stim2(:,stim_idx2)) height(d2.channels)]);
for ii = 1:size(d2.params,2)
    prm = d2.params(:,ii);
   	[~, pred(:,:,ii)] = d2.objFunction(prm, [], stim2(:,stim_idx2), srate);      
end

% Determine time index over which to compute summary statistics
t_idx = d2.t>0.05 & d2.t<1.0;

% Concatenate data and prediction (in order to make sure probabilistic
% assignment of electrodes to areas is done the same way for both).
d = d2.data(t_idx,stim_idx,:);
d = cat(2,d,pred(t_idx, :,:));

% Compute sum across stim_on window
m_conc = squeeze(sum(d,1)); 
[~, channels, group_prob] = groupElecsByVisualArea(d2.channels, 'probabilisticresample', {'V1'});   
[m_conc, se_conc] = averageWithinArea(m_conc, group_prob, [], 10000);

% 1. Contrast response function
subplot('position', posc1); hold on

% Plot linear prediction 
h0 = line([0 x(end)], [0 1], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Plot data
nStim = length(stim_idx);
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'data', 1)

% Plot prediction
m = m_conc(nStim+1:end);
se = se_conc(nStim+1:end,:);
tde_plotPoints(m, se, x2, 'model', 1)

% Format axes
ticklabelsX = num2str(x); ticklabelsX(2:4,:) = ' ';
set(gca, 'xlim', [0 105], 'xtick', x, 'xticklabel', ticklabelsX);
xlabel('Contrast (%)'); ylabel('Summed broadband timecourse (0-1s)'); 
title('Contrast response function', 'fontsize', 18); 

% 2. Peak latency
subplot('position', posc2); hold on

% compute peak across trial window
[~,I] = max(d,[],1);
T = d2.t(t_idx);
[m_conc, se_conc] = averageWithinArea(squeeze(T(I)), group_prob,  [], 10000);

m_conc = m_conc * 1000; % convert to ms
se_conc = se_conc * 1000;

% Plot data
nStim = length(stim_idx);
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'data', 0)

% Plot preduction
m = m_conc(nStim+1:end);
se = se_conc(nStim+1:end,:);
tde_plotPoints(m, se, x2, 'model', 0)

% Format axes
l = get(gca, 'YLim'); ylim([50 l(2)]);
ticklabelsX = num2str(x); ticklabelsX(2:4,:) = ' ';
set(gca, 'xlim', [0 105], 'xtick', x, 'xticklabel', ticklabelsX);

xlabel('Contrast (%)'); ylabel('Peak latency (ms)'); 
title('Latency shift', 'fontsize', 18); 

% 3. Ratio of sustained to  transient
subplot('position', posc3); hold on
T = d2.t(t_idx);

% smooth time courses to get better estimates of max and offset response levels
for ii = 1:size(d,2)
    for jj = 1:size(d,3)
        d(:,ii,jj) = smooth(d(:,ii,jj),40);
    end
end

% DEBUG
%figure; plot(T,squeeze(d(:,5,group_prob>0)));

[M] = max(d,[],1); % value at peak
t_off = (T == 0.5); % offset timepoint
O = d(t_off,:,:); % value at offset
R = squeeze(M./O); % divide value at offset with value of peak

[m_conc, se_conc] = averageWithinArea(R, group_prob,  [], 10000);

% Plot data
nStim = length(stim_idx);
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotSummaryStats(m, se, x, 'data', 0)

% Plot preduction
m = m_conc(nStim+1:end);
se = se_conc(nStim+1:end,:);
tde_plotSummaryStats(m, se, x2, 'model', 0)

% format axes
ticklabelsX = num2str(x); ticklabelsX(2:4,:) = ' ';
set(gca, 'xlim', [0 105], 'xtick', x, 'xticklabel', ticklabelsX);

xlabel('Contrast (%)'); ylabel('Ratio offset - peak'); 
title('Transient-sustained ratio', 'fontsize', 18); 

set(findall(gcf,'-property','FontSize'),'FontSize',20)
