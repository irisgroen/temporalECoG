% tde_mkFigure3

% Load data and fits
modelfun = @DN;
xvalmode = 0;
datatype = 'individualelecs';
[D] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% Select electrodes and compute averages 
[~, ~, group_prob] = groupElecsByVisualArea(D.channels, 'probabilisticresample', {'V1'});   
[data, data_se] = averageWithinArea(D.data, group_prob, @mean, 10000);
[pred, pred_se] = averageWithinArea(D.pred, group_prob, @mean, 10000);
t = D.t;
stim = D.stim;
stim_info = D.stim_info;

% Prepare figure
figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

% Subplot positions: % [left bottom width height]
posa = [0.05 0.55 0.4 0.3];
posb = [0.55 0.5 0.4 0.45];
posc = [0.05 0.05 0.9 0.35];

%% Panel A: example of compressive temporal summation 

% Select two conditions to plot
conditionsOfInterest = {'TWOPULSE-2', 'TWOPULSE-6'};
timepointsOfInterest = [-0.10 1];

stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

s = stim(t_idx,stim_idx);
d = data(t_idx,stim_idx,1);
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
set(gca, 'xticklabel', stim_info.ISI(stim_idx)*1000); box off
xlabel('Stimulus interval (ms)'); ylabel('Neural response'); 
legend({'Stimulus', 'Neural response'}, 'location', 'northwest');
legend('boxoff')

%% Panel B: recovery with adaptation: individual electrodes

% Find stimulus index
stim_idx = find(contains(stim_info.name, {'ONEPULSE-5', 'TWOPULSE'}));
x = stim_info.ISI(stim_idx)*1000; 

% Compute recovery per electrode
srate = D.channels.sampling_frequency(1);
[m] = tde_computeISIrecovery(D.data,D.t,D.stim_info,srate, [], [], 'max');

% Generate new stimuli based on params
nStim = 50;
[stim2, stim_info2] = tde_simulateNewStimuli(t,nStim);
stim_idx2 = find(contains(stim_info2.name, 'TWOPULSE'));
x2 = stim_info2.ISI(stim_idx2)*1000; 

% Predict model responses for new stimuli
pred2 = nan([size(stim2) height(D.channels)]);
for ii = 1:size(D.params,2)
    prm = D.params(:,ii);
   	[~, pred2(:,:,ii)] = D.objFunction(prm, [], stim2, srate);      
end

% Compute recovery per electrode
[m2] = tde_computeISIrecovery_sim(pred2,t,stim_info2,srate);

% Concatenate data and prediction (in order to make sure probabilistic
% assignment of electrodes to areas is done the same way for both).
m = cat(1,m,m2);

% Compute average parameter values within groups
[m_conc, se_conc] = averageWithinArea(m, group_prob, [], 1000);

% Plot
subplot('position', posb); hold on

% Plot linear prediction 
h0 = line([x(1) x(end)], [1 1], 'LineStyle', ':', 'LineWidth', 2, 'color', [0 0 0]);
set(gca, 'xtick', x, 'xticklabelrotation', 45);

% Plot data
nStim = length(stim_idx);
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'errbar', 0,[],50);

% Plot prediction
m = m_conc(nStim+1:end);
se = se_conc(nStim+1:end,:);
tde_plotPoints(m, se, x2, 'ci', 0,[],50,'r');

% Format axes
ylim([0 1.2]);
xlim([-20 x(end)+20]);
xlabel('Stimulus interval (ms)'); ylabel('Ratio second stimulus / first stimulus'); 
legend({'Linear prediction', 'Neural response', 'DN model prediction'}, 'location', 'southeast');
legend('boxoff')

set(findall(gcf,'-property','FontSize'),'FontSize',24)

%% Panel C: data and fits
conditionsOfInterest = {'ONEPULSE-5', 'TWOPULSE'};
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

maxresp = max(d(:,1)); % scale stimulus to max of first condition

subplot('position', posc); hold on
hs = plot(s(:)*maxresp, 'color', [0.7 0.7 0.7], 'linewidth', 1);
hcid = ciplot(d_se_l(:), d_se_u(:), [], 'k', 0.25);
hd = plot(d(:), 'k-', 'linewidth', 2);
hcip = ciplot(p_se_l(:), p_se_u(:), [], 'r', 0.25);
hp = plot(p(:), 'r-', 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(hcid,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(get(get(hcip,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xtick',1:size(d,1):size(d,2)*size(d,1), 'xticklabel', []);
box off;  
ylim([-2 30]);
xlim([-20 length(s(:)) + 20]);

xlabel('Time'); ylabel('Change in broadband (x-fold)'); 
legend({'Neural response', 'DN model prediction'}, 'location', 'northwest');
legend('boxoff');
