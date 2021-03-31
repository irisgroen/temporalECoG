% tde_mkFigure 9

% Load data and fits
modelfun = @DN;%@LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datastr = 'bads';
datatype = 'individualelecs';
D = tde_loadDataForFigure(modelfun, xvalmode, datatype, [], datastr);

% Compute derived parameters
[results] = tde_evaluateModelFit(D);

% Select electrodes and compute averages 
[~, channels, group_prob] = groupElecsByVisualArea(D.channels, 'probabilisticresample');   
[data, data_se] = averageWithinArea(D.data, group_prob, @mean, 1000);
[pred, pred_se] = averageWithinArea(D.pred, group_prob, @mean, 1000);
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

ylabel(sprintf('Response \n (normalized)')); 
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

set(gca, 'xtick',1:size(c1,1):size(c1,2)*size(c1,1), 'xticklabel', stim_info.contrast(stim_idx)*100);
set(gca, 'ytick', [0 1]);
box off,  axis tight

ylim([-0.2 1.5]);
xlim([-20 length(s(:)) + 20]);

xlabel('Contrast (%)'); ylabel(sprintf('Model prediction \n (normalized)')); 
legend(areasOfInterest, 'location', 'northeast');
legend('boxoff');

set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% Panel B: c50 summarized across areas
[m, se] = averageWithinArea(results.derived.params(5,:), group_prob, [], 10000);

subplot('position', posb); hold on
x = 1:height(channels);
tde_plotPoints(m', se, x, 'errbar', 0, [], 40);
set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
ylabel('C50 (%)');
xlabel('Visual area');

% Determine time index over which to compute sum
t_idx = t>0.05 & t<1.0;

% Concatenate data and prediction (in order to make sure probabilistic
% assignment of electrodes to areas is done the same way for both).
d = D.data(t_idx,stim_idx,:);

% Compute sum across stim_on window
sumd = squeeze(max(d,[],1)); 
c50 = [];
for ii = 1:size(sumd,2)
    [~,c50(ii)] = fitNakaRushton(stim_info.contrast(stim_idx)*100,sumd(:,ii));
end

[m, se] = averageWithinArea(c50, group_prob, [], 10000);
tde_plotPoints(m', se, x, 'errbar', 0, [], 40, 'r');
set(findall(gcf,'-property','FontSize'),'FontSize',20)
xlim([0 max(x)+1]); ylim([0 100]);

%% %tmp
[ISIrecover,ts,w] = tde_computeISIrecovery(D.data,D.t,D.stim_info,D.srate,0.5, [], 'max');

[m, se] = averageWithinArea(ts, group_prob, @mean, 1000);

figure;hold on
set(gcf, 'position',  get(0, 'screensize'));
cmap = flipud(brewermap(8,'RdBu'));

for ii = 1:height(channels)
    subplot(2,6,ii);hold on
    p = plot(m(:,:,ii), 'LineWidth', 2);
    title(channels.name(ii));
    axis tight
    set(p, {'color'}, num2cell(cmap,2));

end

[ISIrecover,ts,w] = tde_computeISIrecovery(D.pred,D.t,D.stim_info,D.srate,0.5, [], 'max');


[m, se] = averageWithinArea(ts, group_prob, @mean, 1000);

figure;hold on
set(gcf, 'position',  get(0, 'screensize'));
cmap = flipud(brewermap(8,'RdBu'));

for ii = 1:height(channels)
    subplot(2,6,ii);hold on
    p = plot(m(:,:,ii), 'LineWidth', 2);
    title(channels.name(ii));
    axis tight
    set(p, {'color'}, num2cell(cmap,2));

end
    