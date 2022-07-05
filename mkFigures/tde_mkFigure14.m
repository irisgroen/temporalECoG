% Figure 14 :
% Opposite effects of contrast on response amplitude and timing across
% visual areas.
%
% 2022 Iris Groen

% Load data and fits
modelfun = @DN;
xvalmode = 0;
datatype = 'individualelecs';
datastr = 'fitaverage1000bads';
D = tde_loadDataForFigure(modelfun, xvalmode, datatype, [], datastr);

% Compute derived parameters
[results] = tde_evaluateModelFit(D);

% Select electrodes and compute averages 
if D.options.fitaverage
    % this means electrodes were averaged within area and fitted multiple times
    [data, data_se, channels] = averageMultipleFits(D.data, D.channels, @mean);
    [pred, pred_se] = averageMultipleFits(D.pred, D.channels, @mean);
else
    % this means we have just one fit to all individual electrodes
    numboots = 1000;
    [~, channels, group_prob] = groupElecsByVisualArea(D.channels, 'probabilisticresample');   
    [data, data_se] = averageWithinArea(D.data, group_prob, @mean, numboots);
    [pred, pred_se] = averageWithinArea(D.pred, group_prob, @mean, numboots);
end
t = D.t;
stim_ts = D.stim;
stim_info = D.stim_info;

%% Panels A-E (here figure 1)
% Subplot positions: % [left bottom width height]

posa1 = [0.05  0.75 0.9 0.22];
posa2 = [0.05  0.50 0.9 0.22];
posb =  [0.05  0.1 0.25 0.3];
posc =  [0.375 0.1 0.25 0.3];
posd =  [0.7   0.1 0.25 0.3];

figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

% Panel A: time courses comparison V1 and V3b

conditionsOfInterest = {'CRF'};
timepointsOfInterest = [-0.1 1.2];
areasOfInterest      = {'V1', 'V3b'};

stim_idx  = contains(stim_info.name, conditionsOfInterest);
t_idx     = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);
chan_idx1 = find(strcmp(channels.name, areasOfInterest{1}));
chan_idx2 = find(strcmp(channels.name, areasOfInterest{2}));

% Top: Data

subplot('position', posa1); hold on

s = stim_ts(t_idx,stim_idx);
c1 = data(t_idx,stim_idx,chan_idx1);
c2 = data(t_idx,stim_idx,chan_idx2);
c1 = c1/max(c1(:,end));
c2 = c2/max(c2(:,end));

hs = plot(s(:), 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(c1(:), 'k-', 'linewidth', 2);
hp = plot(c2(:), 'Color', ones(1,3)*0.5, 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xtick',1:size(c1,1):size(c1,2)*size(c1,1), 'xticklabel',[]);
set(gca, 'ytick', [0 1]);
box off,  axis tight

ylim([-0.2 1.5]);
xlim([-20 length(s(:)) + 20]);

ylabel(sprintf('Neural response \n (normalized)')); 
legend(areasOfInterest, 'location', 'northeast');
legend('boxoff');

% Bottom: Model

subplot('position', posa2); hold on

s = stim_ts(t_idx,stim_idx);
c1 = pred(t_idx,stim_idx,chan_idx1);
c2 = pred(t_idx,stim_idx,chan_idx2);
c1 = c1/max(c1(:,end));
c2 = c2/max(c2(:,end));

hs = plot(s(:), 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(c1(:), 'r-', 'linewidth', 2);
hp = plot(c2(:), 'Color', [1 0.5 0.5], 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xtick',1:size(c1,1):size(c1,2)*size(c1,1), 'xticklabel', stim_info.contrast(stim_idx)*100, 'xticklabelrotation', 45);
set(gca, 'ytick', [0 1]);
box off,  axis tight

ylim([-0.2 1.5]);
xlim([-20 length(s(:)) + 20]);

xlabel('Contrast (%)'); ylabel(sprintf('Model prediction \n (normalized)')); 
legend(areasOfInterest, 'location', 'northeast');
legend('boxoff');

set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% Panel B-D: Amplitude differences

% Superimposed normalized time courses

% Select data to plot
conditionsOfInterest = {'CRF-1','CRF-2', 'CRF-3', 'CRF-4','CRF-5'}; 
areasOfInterest= {'V2'};
timepointsOfInterest = [-0.05 0.5];
smoothFactor = 5;

stim_idx = contains(stim_info.name, conditionsOfInterest);
chan_idx = strcmp(channels.name, areasOfInterest);
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

m = data(t_idx,stim_idx,chan_idx);

for ii = 1:size(m,2), m(:,ii) = smooth(m(:,ii),smoothFactor); end

% Make legend
x = stim_info.contrast(stim_idx)*100; %  in percent
l = [];
for ii = length(x):-1:1
    if ii == 1,  l{length(x)+1-ii} = sprintf('%-5.2f', x(ii)); end
    if ii == 2,  l{length(x)+1-ii} = sprintf('%-5.1f', x(ii)); end
    if ii > 2,   l{length(x)+1-ii} = sprintf('%-5.0f', x(ii)); end
end

subplot('position', posb); cla; hold on

colors = gray(length(conditionsOfInterest)+1); colors = colors(1:end-1,:);
ecog_plotMultipleTimeCourses(t(t_idx), fliplr(m), [], colors);
set(gca, 'xlim', timepointsOfInterest);
set(gca, 'ylim', [-0.5 16]);
xlabel('Time (s)'); ylabel('Neural response');

legend(l);
legend('boxoff')

areasOfInterest = {'TO'};
chan_idx = strcmp(channels.name, areasOfInterest);
m = data(t_idx,stim_idx,chan_idx);

for ii = 1:size(m,2), m(:,ii) = smooth(m(:,ii),smoothFactor); end

subplot('position', posc);cla; hold on

ecog_plotMultipleTimeCourses(t(t_idx), fliplr(m), [], colors, [], []);
set(gca, 'xlim', timepointsOfInterest);
set(gca, 'ylim', [-0.06 2.2]);
xlabel('Time (s)'); ylabel('Neural response');

set(findall(gcf,'-property','FontSize'),'FontSize',24)

%% c50

% Data

% Determine time index over which to compute max
t_idx = t>0.05 & t<1.0;

% Get data
d = D.data(t_idx,stim_idx,:);
p = D.pred(t_idx,stim_idx,:);

% Compute max across stim_on window
maxd = squeeze(max(d,[],1));
maxp = squeeze(max(p,[],1));

c50_data  = nan(1,height(D.channels));
c50_model = nan(1,height(D.channels));
fprintf('[%s] Fitting Naka Rushton function to each dataset...\n',mfilename);
for ii = 1:height(D.channels)
    [~,c50_data(ii)] = fitNakaRushton(stim_info.contrast(stim_idx)*100,maxd(:,ii),0);
    [~,c50_model(ii)] = fitNakaRushton(stim_info.contrast(stim_idx)*100,maxp(:,ii),0);
end

if D.options.fitaverage
    [m_data, se_data, channels] = averageMultipleFits(c50_data, D.channels, @mean);
	[m_model, se_model] = averageMultipleFits(c50_model, D.channels, @mean);
else
    [m_data, se_data] = averageWithinArea(c50_data, group_prob, @median, numboots);
    [m_model, se_model] = averageWithinArea(c50_model, group_prob, @median, numboots);
end

% Plot
subplot('position', posd); cla; hold on
x = 1:height(channels);
tde_plotPoints(m_model', se_model, x, 'ci', 0, [], 40, 'r');
tde_plotPoints(m_data', se_data, x, 'errbar', 0, [], 40, 'k');

set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
ylabel('C50 (% contrast)');
xlabel('Visual area');

xlim([0 max(x)+1]); ylim([0 50]);
legend('model', 'data', 'location','northwest');legend boxoff

set(findall(gcf,'-property','FontSize'),'FontSize',24)

%% %% Panel E-F Temporal dynamics differences

figure(2); clf
set(gcf, 'position',  get(0, 'screensize'));

% Select data to plot
conditionsOfInterest = {'CRF-1','CRF-2', 'CRF-3', 'CRF-4','CRF-5'}; 
areasOfInterest= {'V2'};
timepointsOfInterest = [-0.05 0.5];
smoothFactor = 5;

stim_idx = contains(stim_info.name, conditionsOfInterest);
chan_idx = strcmp(channels.name, areasOfInterest);
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

m = data(t_idx,stim_idx,chan_idx);
m = m./max(m);

for ii = 1:size(m,2), m(:,ii) = smooth(m(:,ii),smoothFactor); end

% Make legend
x = stim_info.contrast(stim_idx)*100; %  in percent
l = [];
for ii = length(x):-1:1
    if ii == 1,  l{length(x)+1-ii} = sprintf('%-5.2f', x(ii)); end
    if ii == 2,  l{length(x)+1-ii} = sprintf('%-5.1f', x(ii)); end
    if ii > 2,   l{length(x)+1-ii} = sprintf('%-5.0f', x(ii)); end
end

subplot('position', posb); cla; hold on

colors = gray(length(conditionsOfInterest)+1); colors = colors(1:end-1,:);
ecog_plotMultipleTimeCourses(t(t_idx), fliplr(m), [], colors);
set(gca, 'xlim', timepointsOfInterest);
set(gca, 'ylim', [-0.05 1.05]);
xlabel('Time (s)'); ylabel('Neural response (normalized)');

legend(l);
legend('boxoff')

areasOfInterest = {'TO'};
chan_idx = strcmp(channels.name, areasOfInterest);
m = data(t_idx,stim_idx,chan_idx);
m = m./max(m);

for ii = 1:size(m,2), m(:,ii) = smooth(m(:,ii),smoothFactor); end

subplot('position', posc);cla; hold on
ecog_plotMultipleTimeCourses(t(t_idx), fliplr(m), [], colors, [], []);
set(gca, 'xlim', timepointsOfInterest);
set(gca, 'ylim', [-0.05 1.05]);
xlabel('Time (s)'); ylabel('Neural response (normalized)');

%% Range in time-to-peak

% DATA

% Get data
d = D.data(:,stim_idx,:);
d = d./max(d);

% Find time-to-peak
Id = nan([size(d,2) size(d,3)]);
for ii = 1:size(d,2)
    for jj = 1:size(d,3)
        Id(ii,jj) = find(d(:,ii,jj)==1,1);
    end
end
T = repmat(t, [1 size(d,2) size(d,3)]);
Td = T(Id);
cov_Td = range(Td);

% Average across areas
if D.options.fitaverage
    [T_data, Tse_data, channels] = averageMultipleFits(Td, D.channels, @mean);
	[covT_data, covTse_data] = averageMultipleFits(cov_Td, D.channels, @mean);
else
    [T_data, Tse_data] = averageWithinArea(Td, group_prob, @mean, numboots);
    [covT_data, covTse_data] = averageWithinArea(cov_Td, group_prob, @mean, numboots);
end


% MODEL

% Get predictions
p = D.pred(:,stim_idx,:);
p = p./max(p);

% Find time-to-peak
Ip = nan([size(p,2) size(p,3)]);
for ii = 1:size(p,2)
    for jj = 1:size(p,3)
        %Ip(ii,jj) = find(p(:,ii,jj)>0.5,1);
        Ip(ii,jj) = find(p(:,ii,jj)==1,1);
    end
end
T = repmat(t, [1 size(p,2) size(p,3)]);
Tp = T(Ip);
cov_Tp = range(Tp); 

% Average across areas
if D.options.fitaverage
    [T_pred, Tse_pred, channels] = averageMultipleFits(Tp, D.channels, @mean);
	[covT_pred, covTse_pred] = averageMultipleFits(cov_Tp, D.channels, @mean);
else
    [T_pred, Tse_pred] = averageWithinArea(Tp, group_prob, @median, numboots);
    [covT_pred, covTse_pred] = averageWithinArea(cov_Tp, group_prob, @median, numboots);
end

% Plot range per area
subplot('position', posd); cla; hold on

x = 1:height(channels);
tde_plotPoints(covT_pred', covTse_pred, x, 'ci', 0, [], 40, 'r');
tde_plotPoints(covT_data', covTse_data, x, 'errbar', 0, [], 40, 'k');

set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
ylabel('Range in time-to-peak (s)');
xlabel('Visual area');

xlim([0 max(x)+1]); 
ylim([0 0.55]);
legend('model', 'data', 'location','northwest');legend boxoff

set(findall(gcf,'-property','FontSize'),'FontSize',24)


