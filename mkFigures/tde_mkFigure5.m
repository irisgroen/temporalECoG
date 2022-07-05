% Figure 5 : 
% Comparison of contrast and repetition effects in ECoG data and DN model predictions.
%
% 2022 Iris Groen

modelfun = @DN;
xvalmode = 0;
datatype = 'electrodeaverages';
[D] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

t = D.t;
data = D.data(:,:,1); %V1
pred = D.pred(:,:,1); %V1

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
cmap1    = brewermap(nStim+2, 'Greys');
cmap1    = cmap1(3:end,:);
cmap2    = brewermap(nStim+2, 'OrRd');
cmap2    = cmap2(3:end,:);
xmax = 100; % in ms

% Plot
subplot(2,2,1); hold on
p1 = plot(data(t_idx,stim_idx),'LineWidth', 2);
set(p1, {'color'}, num2cell(cmap1,2));

% Format axes
axis tight
set(gca, 'xlim', [0 xmax], 'ylim', [-2 22], 'xticklabel', []);
xlabel('Time (ms)');
ylabel('Change in power (x-fold)');

% Plot
subplot(2,2,3); hold on
p2 = plot(pred(t_idx,stim_idx),'LineWidth', 2);
set(p2, {'color'}, num2cell(cmap2,2));
xlabel('Time (ms)');
ylabel('Change in power (x-fold)');

% Format axes
axis tight
set(gca, 'xlim', [0 xmax], 'ylim', [-2 22]);
  
%% Right column: adaptation
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
subplot(2,2,2); hold on
p1 = plot(data2(:,1:end-1),'LineWidth', 2);
set(p1, {'color'}, num2cell(cmap1,2));

% Format axes
axis tight
set(gca, 'xlim', [0 xmax], 'ylim', [-2 22], 'xticklabel', [], 'yticklabel', []);

% Plot
subplot(2,2,4); hold on
p2 = plot(pred2(:,1:end-1),'LineWidth', 2);
set(p2, {'color'}, num2cell(cmap2,2));

% Format axes
axis tight
set(gca, 'xlim', [0 xmax], 'ylim', [-2 22], 'yticklabel', []);

set(findall(gcf,'-property','FontSize'),'FontSize',24)

