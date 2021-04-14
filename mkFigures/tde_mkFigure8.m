% tde_mkFigure 8

% Load data and fits
modelfun = @DN;%@LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datastr = 'bads';
datatype = 'individualelecs';
D = tde_loadDataForFigure(modelfun, xvalmode, datatype, [], datastr);

% Compute derived parameters
%[results] = tde_evaluateModelFit(D);

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

conditionsOfInterest = {'ONEPULSE-4', 'TWOPULSE'};
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

ylabel('Response (normalized)'); %title('Broadband responses to increasing durations'); 
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

set(gca, 'xtick',1:size(c1,1):size(c1,2)*size(c1,1), 'xticklabel', stim_info.ISI(stim_idx)*1000, 'xticklabelrotation', 45);
set(gca, 'ytick', [0 1]);
box off,  axis tight

ylim([-0.2 1.5]);
xlim([-20 length(s(:)) + 20]);

xlabel('Stimulus duration (ms)'); %ylabel('Response (normalized)'); %title('Broadband responses to increasing durations'); 
legend(areasOfInterest, 'location', 'northeast');
legend('boxoff');

%% Panel B:recovery for different areas superimposed
[ISIrecover_data] = tde_computeISIrecovery(D.data,D.t,D.stim_info,D.srate,0.5, [], 'max');
[ISIrecover_pred] = tde_computeISIrecovery(D.pred,D.t,D.stim_info,D.srate,0.5, [], 'max');

[m_data, se_data] = averageWithinArea(ISIrecover_data, group_prob, [], 10000);
[m_pred, se_pred] = averageWithinArea(ISIrecover_pred, group_prob, [], 10000);

cmap2      = brewermap(height(channels)+2, '*PuBuGn');
cmap2      = cmap2(1:height(channels),:);

x = stim_info.ISI(stim_idx)*1000; % in ms

chan_ind = 1:5;%height(channels);
subplot('position', posb); hold on

% Plot linear prediction 
h0 = line([x(1) x(end)], [1 1], 'LineStyle', ':', 'LineWidth', 2, 'color', [0 0 0]);
set(gca, 'xtick', x, 'xticklabelrotation', 45);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

for ii = chan_ind
    tde_plotPoints(m_data(:,ii), [], x, 'ci', 0, [],[],cmap2(ii,:));
    %tde_plotPoints(m_pred(:,ii), [], x, 'ci', 0, [],[],cmap2(ii,:));
end
ylim([0 1.5]);
xlim([-20 x(end)+20]);

legend(channels.name(chan_ind), 'location', 'southeast');
legend boxoff
xlabel('Stimulus interval (ms)');
ylabel('Recovery ratio');

% subplot('position', posc); hold on
% 
% % Plot linear prediction 
% h0 = line([x(1) x(end)], [1 1], 'LineStyle', ':', 'LineWidth', 2, 'color', [0 0 0]);
% set(gca, 'xtick', x, 'xticklabelrotation', 45);
% set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
% for ii = chan_ind
%     tde_plotPoints(m_pred(:,ii), [], x, 'ci', 0, [],[],cmap2(ii,:));
%     %tde_plotPoints(m_pred(:,ii), squeeze(se_pred(:,ii,:)), x, 'ci', 0, [],[],cmap2(ii,:));
% end
% ylim([0 1.5]);
% xlim([-20 x(end)+20]);
% xlabel('Stimulus interval (ms)');
% ylabel('Recovery ratio');

%% Panel D: tisi summarized across areas
% [m, se] = averageWithinArea(results.derived.params(4,:), group_prob, [], 10000);
% 
% subplot('position', posd); hold on
% x = 1:height(channels);
% tde_plotPoints(m', se, x, 'errbar', 0, [], 40, 'r');
% xlim([0 max(x)+1]); ylim([0 0.8]);
% set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
% set(gca, 'ytick', 0:0.2:0.8);
% ylabel('Time to recover 95% (s)');
% xlabel('Visual area');

% Generate new stimuli based on params
t = 0/D.srate:1/D.srate:5;
w = 0.5;

% ISI
pulse1_toff = 0.133;
ISIs = 0:(1/D.srate)*5:max(t)-pulse1_toff*2-w;
stimISI = zeros(length(t),length(ISIs)+1);
% generate first pulse
pulse1_on = t>0 & t<=pulse1_toff;

% generate a onepulse condition 
stimISI(pulse1_on,:) = 1;

% generate second pulse
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
[m2, ts2] = tde_computeISIrecovery(pred2,t,stim_info2,D.srate,w,[], 'max');
m2 = round(m2,3);
tISI = [];
for ii = 1:height(D.channels)
    inx = find(m2(:,ii)>0.85,1);
    tISI(ii) = ISIs(inx);
end    

[m, se] = averageWithinArea(tISI, group_prob, [], 10000);

subplot('position', posd); hold on
x = 1:height(channels);
tde_plotPoints(m', se, x, 'errbar', 0, [], 40, 'r');
xlim([0 max(x)+1]); ylim([0 0.6]);
set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
set(gca, 'ytick', 0:0.2:0.6);
ylabel('Time to recover 85% (s)');
xlabel('Visual area');

[m, se] = averageWithinArea(m2, group_prob, [], 10000);
%x = 0:1/D.srate:(size(m,1)-1)*(1/D.srate);
x = ISIs*1000; % in ms

subplot('position', posc); hold on

% % Plot linear prediction 
h0 = line([x(1) x(end)], [1 1], 'LineStyle', ':', 'LineWidth', 2, 'color', [0 0 0]);
set(gca, 'xtick', x, 'xticklabelrotation', 45);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
 
for ii = chan_ind
    %tde_plotPoints(m(:,ii), squeeze(se(:,ii,:)),x, 'ci', 0, [],[],cmap2(ii,:));
    tde_plotPoints(m(:,ii),[],x, 'ci', 0, [],[],cmap2(ii,:));
end
%legend(channels.name); legend boxoff
%legend(channels.name(chan_ind,:));

ylim([0 1.5]);
xlim([-20 600]);
set(gca, 'xtick', stim_info.ISI(stim_idx)*1000); % in ms
xlabel('Stimulus interval (ms)');
ylabel('Recovery ratio');

set(findall(gcf,'-property','FontSize'),'FontSize',20)


% %% %tmp
% [ISIrecover,ts,w] = tde_computeISIrecovery(D.data,D.t,D.stim_info,D.srate,0.5, [], 'max');
% 
% [m, se] = averageWithinArea(ts, group_prob, @mean, 1000);
% 
% figure;hold on
% set(gcf, 'position',  get(0, 'screensize'));
% cmap = flipud(brewermap(8,'RdBu'));
% 
% for ii = 1:height(channels)
%     subplot(2,6,ii);hold on
%     p = plot(m(:,:,ii), 'LineWidth', 2);
%     title(channels.name(ii));
%     axis tight
%     set(p, {'color'}, num2cell(cmap,2));
% 
% end
% 
% [ISIrecover,ts,w] = tde_computeISIrecovery(D.pred,D.t,D.stim_info,D.srate,0.5, [], 'max');
% 
% 
% [m, se] = averageWithinArea(ts, group_prob, @mean, 1000);
% 
% figure;hold on
% set(gcf, 'position',  get(0, 'screensize'));
% cmap = flipud(brewermap(8,'RdBu'));
% 
% for ii = 1:height(channels)
%     subplot(2,6,ii);hold on
%     p = plot(m(:,:,ii), 'LineWidth', 2);
%     title(channels.name(ii));
%     axis tight
%     set(p, {'color'}, num2cell(cmap,2));
% 
% end
    