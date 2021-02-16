% tde_mkFigure3

% Load data and fits

% electrode-averaged data and DN model fits
modelfun = @LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datatype = 'electrodeaverages';
[data1] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% individual electrodes and DN model fits
datatype = 'individualelecs';
[data2] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% subplot positions: % [left bottom width height]

% posa = [0.1 0.55 0.35 0.3];
% posb = [0.1 0.1 0.8 0.25];
% posc = [0.55 0.45 0.35 0.45];

posb = [0.05 0.1 0.9 0.25];
posa = [0.05 0.50 0.25 0.3];
posc = [0.365 0.45 0.25 0.45];
posd = [0.70 0.45 0.25 0.45];

figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

%% Panel A: example of compressive temporal summation 

% Select two conditions to plot
conditionsOfInterest = {'TWOPULSE-2', 'TWOPULSE-6'};
timepointsOfInterest = [-0.10 1];

stim_idx = contains(data1.stim_info.name, conditionsOfInterest);
t_idx    = data1.t>timepointsOfInterest(1) & data1.t<=timepointsOfInterest(2);

s = data1.stim(t_idx,stim_idx);
d = data1.data(t_idx,stim_idx,1);
maxresp = max(d(:,1)); % scale stimulus to max of lowest duration

% for plotting purposes, choose which time points to show
nCut = 240;
nDummy = 60;
s1 = s(1:length(s)-nCut,1);
s1 = cat(1,s1, nan(nDummy,1));
d1 = d(1:length(d)-nCut,1);
d1 = cat(1,d1, nan(nDummy,1));
s2 = s(:,2);
d2 = d(:,2);
sconc = [s1' s2'];
dconc = [d1' d2'];
subplot('position', posa); hold on

plot(sconc, 'color', [0.5 0.5 0.5], 'LineWidth', 1);
plot((dconc./maxresp), 'k', 'lineWidth', 2);
ylim([-0.1 1.4]);
xlim([-20 length(sconc) + 20]);
set(gca, 'xtick',[1 size(s1,1)+1], 'ytick', []);
set(gca, 'xticklabel', data1.stim_info.ISI(stim_idx)*1000); box off
xlabel('Interval (ms)'); ylabel('Neural response'); title('Adaptation to repeated stimulus', 'fontsize', 20); 
legend({'Stimulus', 'Neural data'}, 'location', 'northwest', 'fontsize', 18);
legend('boxoff')

%% Panel B: data and fits
conditionsOfInterest = {'ONEPULSE-5', 'TWOPULSE'};
timepointsOfInterest = [-0.1 1];

stim_idx = contains(data1.stim_info.name, conditionsOfInterest);
t_idx    = data1.t>timepointsOfInterest(1) & data1.t<=timepointsOfInterest(2);

s = data1.stim(t_idx,stim_idx);
d = data1.data(t_idx,stim_idx,1);
p = data1.pred(t_idx,stim_idx,1);
maxresp = max(d(:,1)); % scale stimulus to max of first condition

subplot('position', posb); hold on
hs = plot(s(:)*maxresp, 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(d(:), 'k-', 'linewidth', 2);
hp = plot(p(:), 'r-', 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xtick',1:size(d,1):size(d,2)*size(d,1), 'xticklabel', data2.stim_info.ISI(stim_idx)*1000);
box off;  
ylim([-2 30]);
xlim([-20 length(s(:)) + 20]);

xlabel('Stimulus interval (ms)'); ylabel('Change in broadband (x-fold)'); title('Broadband responses to increasing intervals', 'fontsize', 20); 
legend({'Neural data', 'DN model prediction'}, 'location', 'northwest', 'fontsize', 18);
legend('boxoff');

%% Panel C: recovery with adaptation: timecourses for average

% Compute recovery for average across electrodes
srate = data1.channels.sampling_frequency(1);
[~, ts, w] = tde_computeISIrecovery(data1.data(:,:,1),data1.t,data1.stim_info,srate, [], [], 'sum');

t0 = data1.t(find(data1.t>0,1));
x = t0:(1/srate):w;
x = floor(x*1000); % in ms

subplot('position', posc); hold on
plot(x,ts(:,1), 'k-', 'LineWidth', 4)
colors = flipud(gray(size(ts,2)));
for ii = size(ts,2):-1:2
    plot(x,ts(:,ii), 'Color', colors(ii,:), 'LineWidth', 2);
end

% [~, tspred] = tde_computeISIrecovery(data1.pred(:,:,1),data1.t,data1.stim_info,srate, [], [], 'max');
% subplot('position', posd); hold on
% plot(x,tspred(:,1), 'r:', 'LineWidth', 4)
% colors(:,[2 3]) = 0;
% for ii = size(ts,2):-1:2
%     plot(x,tspred(:,ii), 'LineStyle', ':','Color', colors(ii,:), 'LineWidth', 2);
% end

% make legend
l{1} = 'First stimulus';
stim_idx = find(contains(data2.stim_info.name, {'ONEPULSE-5', 'TWOPULSE'}));
ISIs = data2.stim_info.ISI(stim_idx)*1000; 
for ii = 1:size(ts,2)-1
    l{size(ts,2)+1-ii} = sprintf('ISI %d', ISIs(ii));
end
legend(l, 'fontsize', 18);
legend('boxoff')

% format axes
axis tight
ylim([-2 20]);
xlim([-1 x(end)+1]);
set(gca, 'XTick', [0 x(end)],  'YTick', [0 5 10 15 20]);
xlabel('Time from stim onset (ms)'); ylabel('Response magnitude'); title('Effect of adaptation', 'fontsize', 20); 

%% Panel D: recovery with adaptation: individual electrodes

% Find stimulus index
stim_idx = find(contains(data2.stim_info.name, {'ONEPULSE-5', 'TWOPULSE'}));
x = data2.stim_info.ISI(stim_idx)*1000; 

% Compute recovery per electrode
srate = data2.channels.sampling_frequency(1);
[m] = tde_computeISIrecovery(data2.data,data2.t,data2.stim_info, srate, [], [], 'sum');

% Generate new stimuli based on params
nStim = 50;
[stim2, stim_info2] = tde_simulateNewStimuli(data2.t,nStim);
stim_idx2 = find(contains(stim_info2.name, 'TWOPULSE'));
x2 = stim_info2.ISI(stim_idx2)*1000; 

% Predict model responses for new stimuli
srate = data2.channels.sampling_frequency(1);
pred = nan([size(stim2) height(data2.channels)]);
for ii = 1:size(data2.params,2)
    prm = data2.params(:,ii);
   	[~, pred(:,:,ii)] = data2.objFunction(prm, [], stim2, srate);      
end

% Compute recovery per electrode
[m2] = tde_computeISIrecovery_sim(pred,data2.t,stim_info2,srate);

% Concatenate data and prediction (in order to make sure probabilistic
% assignment of electrodes to areas is done the same way for both).
m = cat(1,m,m2);

% Compute average parameter values within groups
[~, channels, group_prob] = groupElecsByVisualArea(data2.channels, 'probabilisticresample', {'V1'});   
[m_conc, se_conc] = averageWithinArea(m, group_prob, [], 10000);

% Plot
subplot('position', posd); hold on

% Plot linear prediction 
h0 = line([x(1) x(end)], [1 1], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);

% Plot data
nStim = length(stim_idx);
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'data', 0)

% Plot prediction
m = m_conc(nStim+1:end);
se = se_conc(nStim+1:end,:);
tde_plotPoints(m, se, x2, 'model', 0)

% Format axes
ylim([0 1.2]);
xlim([-20 x(end)+20]);
xlabel('Stimulus interval (ms)'); ylabel('Ratio second stimulus / first stimulus'); title('Recovery from adaptation', 'fontsize', 20); 
legend({'Linear prediction', 'Neural data', 'DN model prediction'}, 'location', 'southeast', 'fontsize', 18);
legend('boxoff')

set(findall(gcf,'-property','FontSize'),'FontSize',20)
