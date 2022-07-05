% Figure 6 : 
% How the DN model predicts similar effects of low-contrast and repetition.
%
% 2022 Iris Groen

modelfun = @DN;
xvalmode = 0;
datatype = 'electrodeaverages';
[D] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

t = D.t;
param = D.params(:,1); %V1
data = D.data(:,:,1); %V1
pred = D.pred(:,:,1); %V1
stim = D.stim;
srate = D.srate;

[~,~, numrsp,demrsp] = DN_modifiedIRF(param, data, stim, srate);

% Prepare figure
figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

%% Left column: contrast

conditionsOfInterest = {'CRF'};
timepointsOfInterest = [0 0.5];

stim_idx = find(contains(D.stim_info.name, conditionsOfInterest));
nStim    = length(stim_idx);
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

% Get colormaps for plotting
cmap2    = brewermap(nStim+2, 'OrRd');
cmap2    = cmap2(3:end,:);

conditionsOfInterest = {'CRF-5'};
cond_idx = find(contains(D.stim_info.name(stim_idx), conditionsOfInterest));

subplot(2,2,1); cla; hold on
hs = plot(t,stim(:,stim_idx(cond_idx)), 'color', [0.7 0.7 0.7], 'linewidth', 1);
plot(t,numrsp(:,stim_idx(cond_idx)),'LineWidth', 2, 'Color', cmap2(cond_idx,:));
plot(t,demrsp(:,stim_idx(cond_idx)),'LineWidth', 2, 'Color', cmap2(cond_idx,:), 'LineStyle', '--');
%plot(t,numrsp(:,stim_idx(cond))./demrsp(:,stim_idx(cond)),'LineWidth', 2, 'Color', cmap2(4,:));
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
legend({'Numerator', 'Denominator'}, 'location', 'northwest'); legend boxoff

set(gca, 'xlim', [-0.1 0.65], 'ylim', [0 1.1]);
ylabel('Model output (a.u.)');

conditionsOfInterest = {'CRF-1'};
cond_idx = find(contains(D.stim_info.name(stim_idx), conditionsOfInterest));

subplot(2,2,3); cla; hold on
hs = plot(t,stim(:,stim_idx(cond_idx)), 'color', [0.7 0.7 0.7], 'linewidth', 1);
plot(t,numrsp(:,stim_idx(cond_idx)),'LineWidth', 2, 'Color', cmap2(cond_idx,:));
plot(t,demrsp(:,stim_idx(cond_idx)),'LineWidth', 2, 'Color', cmap2(cond_idx,:), 'LineStyle', '--');
%plot(t,numrsp(:,stim_idx(cond))./demrsp(:,stim_idx(cond)),'LineWidth', 2, 'Color', cmap2(4,:));
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xlim', [-0.1 0.65], 'ylim', [0 0.2]);
xlabel('Time (ms)');
ylabel('Model output (a.u.)');

  
%% Right column: adaptation
w = 0.3;
[~, data2] = tde_computeISIrecovery(data,D.t,D.stim_info,D.srate,w);
[~, pred2] = tde_computeISIrecovery(pred,D.t,D.stim_info,D.srate,w);
Ts = {'ONEPULSE-5', 'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6','first pulse'};
nStim = length(Ts);

% Get colormaps for plotting
cmap2    = brewermap(nStim+2, 'OrRd');
cmap2    = cmap2(3:end,:);

conditionsOfInterest = {'TWOPULSE-5'};
stim_idx = find(contains(D.stim_info.name, conditionsOfInterest));
cond_idx = find(contains(Ts, conditionsOfInterest));
 
% Plot
subplot(2,2,2); cla; hold on
hs = plot(t,stim(:,stim_idx), 'color', [0.7 0.7 0.7], 'linewidth', 1);
plot(t,numrsp(:,stim_idx),'LineWidth', 2, 'Color', cmap2(cond_idx,:));
plot(t,demrsp(:,stim_idx),'LineWidth', 2, 'Color', cmap2(cond_idx,:), 'LineStyle', '--');
%plot(t,numrsp(:,stim_idx)./demrsp(:,stim_idx),'LineWidth', 2, 'Color', cmap2(4,:));
%legend({'numerator', 'denominator'}, 'location', 'northeast'); legend boxoff
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xlim', [-0.1 0.65], 'ylim', [0 1]);

conditionsOfInterest = {'TWOPULSE-2'};
stim_idx = find(contains(D.stim_info.name, conditionsOfInterest));
cond_idx = find(contains(Ts, conditionsOfInterest));

% Plot
subplot(2,2,4); cla; hold on
hs = plot(t,stim(:,stim_idx), 'color', [0.7 0.7 0.7], 'linewidth', 1);
plot(t,numrsp(:,stim_idx),'LineWidth', 2, 'Color', cmap2(cond_idx,:));
plot(t,demrsp(:,stim_idx),'LineWidth', 2, 'Color', cmap2(cond_idx,:), 'LineStyle', '--');
%plot(t,numrsp(:,stim_idx)./demrsp(:,stim_idx),'LineWidth', 2, 'Color', cmap2(4,:));
%legend({'numerator', 'denominator'}, 'location', 'northeast'); legend boxoff
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xlim', [-0.1 0.65], 'ylim', [0 1]);
xlabel('Time (ms)');

set(findall(gcf,'-property','FontSize'),'FontSize',24)

