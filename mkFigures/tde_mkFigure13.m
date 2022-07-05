% Figure 13 :
% Temporal adaptation windows do not differ systematically between visual
% areas.
% 
% 2022 Iris Groen

% Load data and fits
modelfun = @DN;
xvalmode = 0;
datatype = 'individualelecs';
datastr = 'fitaverage1000bads';
D = tde_loadDataForFigure(modelfun, xvalmode, datatype, [], datastr);

% Note: last panel takes a long time to plot!!

%% Panel A: Example time courses for ISI stimuli

% Select electrodes and compute averages 
if D.options.fitaverage
    % this means electrodes were averaged within area and fitted multiple times
    [data, data_se, channels] = averageMultipleFits(D.data, D.channels, @mean);
    [pred, pred_se] = averageMultipleFits(D.pred, D.channels, @mean);
else
    % this means we have just one fit to all individual electrodes 
    numboot = 1000;
    [~, channels, group_prob] = groupElecsByVisualArea(D.channels, 'probabilisticresample');   
    [data, data_se] = averageWithinArea(D.data, group_prob, @mean, numboot);
    [pred, pred_se] = averageWithinArea(D.pred, group_prob, @mean, numboot);
end
t = D.t;
stim_ts = D.stim;
stim_info = D.stim_info;

% Plot specs
% Subplot positions: % [left bottom width height]

posa1 = [0.05  0.75 0.9 0.22];
posa2 = [0.05  0.50 0.9 0.22];
posb =  [0.05  0.1 0.25 0.3];
posc =  [0.375 0.1 0.25 0.3];
posd =  [0.7   0.1 0.25 0.3];

figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

% Panel A: time courses comparison V1 and V3b

conditionsOfInterest = {'ONEPULSE-4', 'TWOPULSE'};
timepointsOfInterest = [-0.1 1.2];
areasOfInterest      = {'V1', 'V3b'};

% Get timecourses
stim_idx  = contains(stim_info.name, conditionsOfInterest);
t_idx     = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);
chan_idx1 = find(strcmp(channels.name, areasOfInterest{1}));
chan_idx2 = find(strcmp(channels.name, areasOfInterest{2}));

% Top: Data

subplot('position', posa1);cla; hold on

s = stim_ts(t_idx,stim_idx);
c1 = data(t_idx,stim_idx,chan_idx1);
c2 = data(t_idx,stim_idx,chan_idx2);
c1 = c1/max(c1(:,end));
c2 = c2/max(c2(:,end));

% Plot
hs = plot(s(:), 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(c1(:), 'k-', 'linewidth', 2);
hp = plot(c2(:), 'Color', ones(1,3)*0.5, 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Format axes
set(gca, 'xtick',1:size(c1,1):size(c1,2)*size(c1,1), 'xticklabel',[]);
set(gca, 'ytick', [0 1]);
box off,  axis tight

ylim([-0.2 1.5]);
xlim([-20 length(s(:)) + 20]);

ylabel(sprintf('Neural response \n (normalized)')); 
legend(areasOfInterest, 'location', 'northeast');
legend('boxoff');

% Bottom: Model

subplot('position', posa2); cla; hold on

% Get timecourses
s = stim_ts(t_idx,stim_idx);
c1 = pred(t_idx,stim_idx,chan_idx1);
c2 = pred(t_idx,stim_idx,chan_idx2);
c1 = c1/max(c1(:,end));
c2 = c2/max(c2(:,end));

% Plot
hs = plot(s(:), 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(c1(:), 'r-', 'linewidth', 2);
hp = plot(c2(:), 'Color', [1 0.5 0.5], 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Format axes
set(gca, 'xtick',1:size(c1,1):size(c1,2)*size(c1,1), 'xticklabel', stim_info.ISI(stim_idx)*1000, 'xticklabelrotation', 45);
set(gca, 'ytick', [0 1]);
box off,  axis tight

ylim([-0.2 1.5]);
xlim([-20 length(s(:)) + 20]);

ylabel(sprintf('Model prediction \n (normalized)')); 
xlabel('Stimulus duration (ms)'); %ylabel('Response (normalized)'); %title('Broadband responses to increasing durations'); 
legend(areasOfInterest, 'location', 'northeast');
legend('boxoff');

%% Panel B:recovery for different areas superimposed: data
[ISIrecover_data] = tde_computeISIrecovery(D.data,D.t,D.stim_info,D.srate,0.5, [], 'max');

if D.options.fitaverage
    [m_data,~, channels] = averageMultipleFits(ISIrecover_data, D.channels, @mean);
else
    [m_data] = averageWithinArea(ISIrecover_data, group_prob, @median, numboot);
end

cmap      = brewermap(height(channels)+2, '*PuBuGn');
cmap      = cmap(1:height(channels),:);

x = stim_info.ISI(stim_idx)*1000; % in ms

chan_ind = height(channels):-1:1;
subplot('position', posb); cla; hold on

% Plot linear prediction 
h0 = line([x(1) x(end)], [1 1], 'LineStyle', ':', 'LineWidth', 2, 'color', [0 0 0]);
set(gca, 'xtick', x, 'xticklabelrotation', 45);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

for ii = chan_ind
    tde_plotPoints(m_data(:,ii), [], x, 'ci', 0, [],[],cmap(ii,:));
end

% Format axes
ylim([0 1.1]);
xlim([-20 x(end)+20]);

legend(channels.name(chan_ind), 'location', 'southeast');
legend boxoff
xlabel('Stimulus interval (ms)');
ylabel('Recovery ratio');

%% Panel C: recovery for multiple areas superimposed: model

% Generate new stimuli based on params
t = 0/D.srate:1/D.srate:2;
w = 0.5;

% ISI
pulse1_toff = 0.133;
ISIs = 0:(1/D.srate)*5:max(t)-pulse1_toff*2-w;
stimISI = zeros(length(t),length(ISIs)+1);

% Generate first pulse
pulse1_on = t>0 & t<=pulse1_toff;

% Generate a ONEPULSE condition 
stimISI(pulse1_on,:) = 1;

% Generate second pulse
for ii = 1:length(ISIs)
    this_ISI = ISIs(ii);%ii/nStim * 0.533;
    pulse2_ton = pulse1_toff + this_ISI;
    pulse2_toff = pulse2_ton + 0.133;
    pulse2_on = t>pulse2_ton & t<=pulse2_toff;
    stimISI(pulse2_on,ii+1) = 1;
end

% Predict model responses for new stimuli
pred2 = nan([size(stimISI) height(D.channels)]);
for ii = 1:size(D.params,2)
    prm = D.params(:,ii);
   	[~, pred2(:,:,ii)] = D.objFunction(prm, [], stimISI, D.srate);      
end

nStim = size(stimISI,2);
name = ['ONEPULSE-4'; repmat({'TWOPULSE'}, nStim-1,1)];
ISI = [0; ISIs'];
duration = ones(nStim,1)*0.133;
stim_info2 = table(name, duration, ISI);

% Compute recovery per electrode
[ISIrecover_pred] = tde_computeISIrecovery(pred2,t,stim_info2,D.srate,w,[], 'max');

if D.options.fitaverage
    [m_pred,~,channels] = averageMultipleFits(ISIrecover_pred, D.channels, @mean);
else
    [m_pred] = averageWithinArea(ISIrecover_pred, group_prob, @median, numboot);
end

cmap      = brewermap(height(channels)+2, '*YlOrRd');
cmap      = cmap(1:height(channels),:);

% Plot
subplot('position', posc); cla; hold on

% % Plot linear prediction 
h0 = line([x(1) x(end)], [1 1], 'LineStyle', ':', 'LineWidth', 2, 'color', [0 0 0]);
set(gca, 'xtick', x, 'xticklabelrotation', 45);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
 
for ii = chan_ind
    tde_plotPoints(m_pred(:,ii),[],ISIs*1000, 'ci', 0, [],[],cmap(ii,:));
end

% Format axes
ylim([0 1.1]);
xlim([-20 x(end)+20]);
set(gca, 'xtick', stim_info.ISI(stim_idx)*1000); % in ms
xlabel('Stimulus interval (ms)');
ylabel('Recovery ratio');
legend(channels.name(chan_ind), 'location', 'southeast');
legend boxoff

%% Panel D: tisi summarized across areas for both data and model

% Set a threshold to measure recovery at (%)
thresh = 0.80; 

x = stim_info.ISI(stim_idx)*1000; % in ms

% DATA
fprintf('[%s] Computing tISI for data...\n',mfilename);

% Fit an exponential to the recovery of the data points to obtain an estimate of
% the recovery over time
formula_to_fit = 'a * x ^ b + c';
sp = [1 1 1]; lb = [0 0 0]; ub = [1 1 1]; f = cell(1,height(D.channels));
for ii = 1:height(D.channels)
	y = ISIrecover_data(:,ii);
	f{ii} = fit(x, y, formula_to_fit, 'StartPoint', sp, 'Lower', lb, 'Upper', ub);
end

% Generate the prediction based on the fitted line
x1 = 0:max(x)/1000:max(x)*4; 
y1 = nan(length(x1), height(D.channels));
for ii = 1:height(D.channels)
    f1 = f{ii};
    y1(:,ii) = f1(x1);
end

% Find the ISI at which the lines crosses the threshold; this is tISI
% derived from the data
tISI_data = nan(1,height(D.channels));
for ii = 1:height(D.channels)
    inx = find(y1(:,ii)>thresh,1);
    if isempty(inx)
        inx = length(x1);
    end
    tISI_data(ii) = x1(inx);
end

% Average across electrodes or fits
if D.options.fitaverage
    [m_data, se_data, channels] = averageMultipleFits(tISI_data, D.channels, @mean);
else
    [m_data, se_data] = averageWithinArea(tISI_data, group_prob, @median, numboot);
end

% MODEL

fprintf('[%s] Computing tISI for model...\n',mfilename);

% Find the ISI at which the predicted recovery from the model crosses the
% threshold; this is tISI derived from the model
ISIrecover_pred = round(ISIrecover_pred,3);
tISI_pred = nan(1,height(D.channels));
for ii = 1:height(D.channels)
    inx = find(ISIrecover_pred(:,ii)>thresh,1);
    if isempty(inx)
        inx = size(ISIrecover_pred,1);
    end
    tISI_pred(ii) = ISIs(inx);
end    

% Average across electrodes or fits
if D.options.fitaverage
    [m_pred, se_pred, channels] = averageMultipleFits(tISI_pred, D.channels, @mean);
else
    [m_pred, se_pred] = averageWithinArea(tISI_pred, group_prob, @median, numboot);
end

% Plot
subplot('position', posd); cla; hold on
x = 1:height(channels);
tde_plotPoints(m_pred', se_pred, x, 'ci', 0, [], 40, 'r');
tde_plotPoints((m_data/1000)',se_data/1000, x, 'errbar', 0, [], 40, 'k');

% Format axes
xlim([0 max(x)+1]); ylim([0 0.8]);
set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
set(gca, 'ytick', 0:0.2:1);
ylabel(sprintf('Time to recover %d percent (s)', thresh*100));
xlabel('Visual area');

legend('model', 'data', 'Location', 'NorthWest')
legend boxoff

set(findall(gcf,'-property','FontSize'),'FontSize',24)


    