% tde_mkFigure5

% contrast response function + second pulse adaptation 
% next to DN model predictions

modelfun = @LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datatype = 'electrodeaverages';
[D] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% Prepare figure
figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

t = D.t;
param = D.params(:,1);
data = D.data(:,:,1);
pred = D.pred(:,:,1);
stim = D.stim;
srate = D.srate;

[~,~, numrsp,demrsp] = LINEAR_RECTF_EXP_NORM_DELAY(param, data, stim, srate);

%% Top row: contrast
clf
conditionsOfInterest = {'CRF'};
timepointsOfInterest = [0 0.5];

stim_idx = find(contains(D.stim_info.name, conditionsOfInterest));
nStim    = length(stim_idx);
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

% Get colormaps for plotting
cmap1    = brewermap(nStim+2, 'Greys');
cmap1    = cmap1(3:end,:);
cmap2    = brewermap(nStim+2, 'OrRd');
cmap2    = cmap2(3:end,:);
xmax = 100;

% Plot
subplot(2,3,1); hold on
p1 = plot(data(t_idx,stim_idx),'LineWidth', 2);
set(p1, {'color'}, num2cell(cmap1,2));

% Format axes
axis tight
set(gca, 'xlim', [0 xmax], 'ylim', [-2 22]);
xlabel('Time (ms)');
ylabel('Change in power (xfold)');

subplot(2,3,2); hold on
p2 = plot(pred(t_idx,stim_idx),'LineWidth', 2);
set(p2, {'color'}, num2cell(cmap2,2));

% Format axes
axis tight
set(gca, 'xlim', [0 xmax], 'ylim', [-2 22]);
xlabel('Time (ms)');
ylabel('Change in power (xfold)');

timepointsOfInterest = [0 1];
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);
cond = 3;

subplot(2,3,3); hold on
hs = plot(t,stim(:,stim_idx(cond)), 'color', [0.7 0.7 0.7], 'linewidth', 1);
plot(t,numrsp(:,stim_idx(cond)),'LineWidth', 2, 'Color', 'b');
plot(t,demrsp(:,stim_idx(cond)),'LineWidth', 2, 'Color', 'c');
%plot(numrsp(t_idx,stim_idx(cond))./demrsp(t_idx,stim_idx(cond)),'LineWidth', 2, 'Color', cmap2(4,:));
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

legend({'Numerator', 'Denominator'}, 'location', 'northeast'); legend boxoff
set(gca, 'xlim', [-0.15 0.6], 'ylim', [0 0.4]);
xlabel('Time (ms)');
ylabel('Model output (a.u.)');
  
%% Bottom row: adaptation
w = 0.3;
[~, data2] = tde_computeISIrecovery(data,D.t,D.stim_info,D.srate,w);
[~, pred2] = tde_computeISIrecovery(pred,D.t,D.stim_info,D.srate,w);
Ts = {'ONEPULSE-5', 'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6','first pulse'};

% Get colormaps for plotting
nStim    = size(pred2,2)-1;
cmap1    = brewermap(nStim+2, 'Greys');
cmap1    = cmap1(3:end,:);
cmap2    = brewermap(nStim+2, 'OrRd');
cmap2    = cmap2(3:end,:);
xmax     = 100;

% Plot
subplot(2,3,4); hold on
p1 = plot(data2(:,1:end-1),'LineWidth', 2);
set(p1, {'color'}, num2cell(cmap1,2));

% Format axes
axis tight
set(gca, 'xlim', [0 xmax], 'ylim', [-2 22]);
xlabel('Time (ms)');
ylabel('Change in power (xfold)');

% Plot
subplot(2,3,5); hold on
p2 = plot(pred2(:,1:end-1),'LineWidth', 2);
set(p2, {'color'}, num2cell(cmap2,2));

% Format axes
axis tight
set(gca, 'xlim', [0 xmax], 'ylim', [-2 22]);
xlabel('Time (ms)');
ylabel('Change in power (xfold)');

conditionsOfInterest = {'TWOPULSE-3'};
stim_idx = find(contains(D.stim_info.name, conditionsOfInterest));
timepointsOfInterest = [0 1];
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

% Plot
subplot(2,3,6); hold on
cond = 4;
hs = plot(t,stim(:,stim_idx), 'color', [0.7 0.7 0.7], 'linewidth', 1);
plot(t,numrsp(:,stim_idx),'LineWidth', 2, 'Color', 'b');
plot(t,demrsp(:,stim_idx),'LineWidth', 2, 'Color', 'c');
%plot(t,numrsp(t_idx,stim_idx(cond))./demrsp(t_idx,stim_idx(cond)),'LineWidth', 2, 'Color', cmap2(4,:));
%legend({'numerator', 'denominator'}, 'location', 'northeast'); legend boxoff
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xlim', [0.1 0.6], 'ylim', [0 1.1]);
xlabel('Time (ms)');
ylabel('Model output (a.u.)');

set(findall(gcf,'-property','FontSize'),'FontSize',24)

