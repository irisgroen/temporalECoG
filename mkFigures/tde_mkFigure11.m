% Figure 11 : 
% Temporal summation window sizes increase in higher visual areas.
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
    % This means electrodes were averaged within area and fitted multiple times
    [data, data_se, channels] = averageMultipleFits(D.data, D.channels, @mean);
    [pred, pred_se] = averageMultipleFits(D.pred, D.channels, @mean);
else
    % This means we have just one fit to each of the individual electrodes 
    numboot = 1000;
    [~, channels, group_prob] = groupElecsByVisualArea(D.channels, 'probabilisticresample');   
    [data, data_se] = averageWithinArea(D.data, group_prob, @mean, numboot);
    [pred, pred_se] = averageWithinArea(D.pred, group_prob, @mean, numboot);
end
t = D.t;
stim_ts = D.stim;
stim_info = D.stim_info;

%% Plot specs
% Subplot positions: % [left bottom width height]

posa1 = [0.05  0.75 0.9 0.22];
posa2 = [0.05  0.50 0.9 0.22];
posb =  [0.05  0.1 0.25 0.3];
posc =  [0.375 0.1 0.25 0.3];
posd =  [0.7   0.1 0.25 0.3];

figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

% Panel A: time courses comparison V1 and V3b

COI = {'ONEPULSE'};
TOI = [-0.1 0.8];
AOI = {'V1', 'V3b'};

stim_idx  = contains(stim_info.name, COI);
t_idx     = t>TOI(1) & t<=TOI(2);
chan_idx1 = find(strcmp(channels.name, AOI{1}));
chan_idx2 = find(strcmp(channels.name, AOI{2}));

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
legend(AOI, 'location', 'northeast');
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

set(gca, 'xtick',1:size(c1,1):size(c1,2)*size(c1,1), 'xticklabel', stim_info.duration(stim_idx)*1000, 'xticklabelrotation', 45);
set(gca, 'ytick', [0 1]);
box off,  axis tight

ylim([-0.2 1.5]);
xlim([-20 length(s(:)) + 20]);

ylabel(sprintf('Model prediction \n (normalized)')); 
xlabel('Stimulus duration (ms)'); 
legend(AOI, 'location', 'northeast');
legend('boxoff');

%% Panel B, C, D: derived parameters summarized across areas

% Data
[derivedPrm] = tde_computeDerivedParamsData(D.data(:,stim_idx,:),t); 

if D.options.fitaverage
    [m_data, se_data, channels] = averageMultipleFits(derivedPrm, D.channels, @mean);
else
    [m_data, se_data] = averageWithinArea(derivedPrm, group_prob, @median, numboot);
end

% Model
if D.options.fitaverage
    [m_model, se_model, channels] = averageMultipleFits(results.derived.params, D.channels, @mean);
else
    [m_model, se_model] = averageWithinArea(results.derived.params, group_prob, @median, numboot);
end

subplot('position', posb); cla; hold on
x = 1:height(channels);
tde_plotPoints(m_model(1,:)', squeeze(se_model(1,:,:)), x, 'ci', 0, [], 40, 'r');
tde_plotPoints(m_data(1,:)', squeeze(se_data(1,:,:)), x, 'errbar', 0, [], 40);
xlim([0 max(x)+1]); ylim([0 0.2]);
set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
set(gca, 'ytick', [0 0.1 0.2]);
ylabel('Time to peak (s)');
xlabel('Visual area');
legend('model', 'data', 'Location', 'NorthWest'); legend boxoff

subplot('position', posc); cla; hold on
x = 1:height(channels);
tde_plotPoints(m_model(3,:)', squeeze(se_model(3,:,:)), x, 'ci', 0, [], 40, 'r');
tde_plotPoints(m_data(3,:)', squeeze(se_data(3,:,:)), x, 'errbar', 0, [], 40);

xlim([0 max(x)+1]); ylim([0 0.2]);
set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
set(gca, 'ytick', [0 0.1 0.2]);
ylabel('Fwhm (s)');
xlabel('Visual area');

subplot('position', posd); cla; hold on
x = 1:height(channels);
tde_plotPoints(m_model(2,:)', squeeze(se_model(2,:,:)), x, 'ci', 0, [], 40, 'r');
tde_plotPoints(m_data(2,:)', squeeze(se_data(2,:,:)), x, 'errbar', 0, [], 40);
xlim([0 max(x)+1]); ylim([0 1]);
set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
set(gca, 'ytick', 0:0.2:1);
ylabel('Ratio sustained/transient');
xlabel('Visual area');

set(findall(gcf,'-property','FontSize'),'FontSize',24)
