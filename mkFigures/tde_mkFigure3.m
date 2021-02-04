% tde_mkFigure3

% Load data and fits

% electrode-averaged data and DN model fits
data1 = load('/Users/iiagroen/surfdrive/BAIR/Papers/TemporalDynamicsECoG/results/LINEAR_RECTF_EXP_NORM_DELAY_xvalmode0_electrodeaverages');

% individual electrodes and DN model fits
data2 = load('/Users/iiagroen/surfdrive/BAIR/Papers/TemporalDynamicsECoG/results/LINEAR_RECTF_EXP_NORM_DELAY_xvalmode0_individualelecs');

% subplot positions: % [left bottom width height]
posa = [0.1 0.55 0.35 0.3];
posb = [0.1 0.1 0.8 0.25];
posc = [0.55 0.45 0.35 0.45];

%% Panel A: example of compressive temporal summation 
%figure(1); clf
figure; hold on;
set(gcf, 'position',  get(0, 'screensize'));

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
ylim([-0.2 1.4]);
xlim([-20 length(sconc) + 20]);
set(gca, 'xtick',[1 size(s1,1)+1]);
set(gca, 'xticklabel', data1.stim_info.ISI(stim_idx)*1000); box off
xlabel('interval (ms)'); ylabel('neural response'); title('Adaptation to repeated stimulus', 'fontsize', 20); 
legend({'stimulus', 'neural data'}, 'location', 'northwest', 'fontsize', 18);
legend('boxoff')

%% Panel B: data and fits
conditionsOfInterest = {'TWOPULSE'};
timepointsOfInterest = [-0.1 1];

stim_idx = contains(data1.stim_info.name, conditionsOfInterest);
t_idx    = data1.t>timepointsOfInterest(1) & data1.t<=timepointsOfInterest(2);

s = data1.stim(t_idx,stim_idx);
d = data1.data(t_idx,stim_idx,1);
p = data1.pred(t_idx,stim_idx,1);

subplot('position', posb); hold on
hs = plot(s(:), 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(d(:), 'k-', 'linewidth', 2);
hp = plot(p(:), 'r-', 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


set(gca, 'xtick',1:size(d,1):size(d,2)*size(d,1), 'xticklabel', data2.stim_info.ISI(stim_idx)*1000);
box off;  
ylim([-2 30]);
xlim([-20 length(s(:)) + 20]);

xlabel('stimulus interval (ms)'); ylabel('x-fold change in broadband'); title('V1 timecourses and fits', 'fontsize', 20); 
legend({'neural data', 'DN model prediction'}, 'location', 'northwest', 'fontsize', 18);
legend('boxoff');

%% Panel C: recovery with adaptation
stim_inx = find(contains(data2.stim_info.name, {'ONEPULSE-5', 'TWOPULSE'}));

stim_on = data2.t>0; % compute sum over entire epoch, so we don't miss the second pulse; 
x = data2.stim_info.ISI(stim_inx)*1000; 

subplot('position', posc); hold on

% Compute recovery per electrode
[m] = tde_computeISIrecovery(data2.data,data2.t,data2.stim_info);
%[m] = tde_computeISIrecovery(data1.data(:,:,1),data1.t,data1.stim_info);

%%
% Compute average parameter values within groups
[~, channels, group_prob] = groupElecsByVisualArea(data2.channels, 'probabilisticresample', {'V1'});   
[m, se] = averageWithinArea(m, group_prob, [], 10000);

% Plot linear prediction 
h0 = line([x(1) x(end)], [100 100], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);

% Plot averages and fit
formula_to_fit = 'a * x ^ b + c';
sp = [1 1 1];
[f] = fitData(x,m,1,formula_to_fit, sp);
[y] = plotData(x,m,se,f,1,[0 0 0],0);  

% Generate new stimuli based on params
nStim = 50;
[stim, stim_info] = tde_simulateNewStimuli(data2.t,nStim);

% Predict model responses for new stimuli
srate = data2.channels.sampling_frequency(1);
pred = nan([size(stim) height(data2.channels)]);
for ii = 1:size(data2.params,2)
    prm = data2.params(:,ii);
   	[~, pred(:,:,ii)] = data2.objFunction(prm, [], stim, srate);      
end

[m2] = tde_computeISIrecovery_sim(pred,data2.t,stim_info);
[m2, se2] = averageWithinArea(m2, group_prob, [], 10000);

% % Plot average model predictions across electrodes
stim_inx = find(contains(stim_info.name, 'TWOPULSE'));
x2 = stim_info.ISI(stim_inx)*1000;
plot(x2,m2, 'r', 'linewidth', 2);
ch = ciplot(se2(:,1), se2(:,2), x2, 'r', 0.25);
set(get(get(ch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
 
% format axes
ylim([0 120]);
xlim([-20 x(end)+20]);
xlabel('stimulus interval (ms)'); ylabel('second stimulus % of first stimulus'); title('Recovery from adaptation in V1', 'fontsize', 20); 
%legend({'linear prediction', 'neural data', 'fitted line', 'DN model prediction'}, 'location', 'southeast', 'fontsize', 18);
legend({'linear prediction', 'neural data', 'DN model prediction'}, 'location', 'southeast', 'fontsize', 18);
legend('boxoff')
%axis square
%axis tight

%% subfunctions

function [f] = fitData(x,m,chan_plot_idx,formula_to_fit, sp, lb, ub)
    if ~exist('sp', 'var'), sp = []; end
    if ~exist('lb', 'var'), lb = []; end
    if ~exist('ub', 'var'), ub = []; end
    for ii = 1:length(chan_plot_idx)
        y = m(:,chan_plot_idx(ii));
        f{ii} = fit(x, y, formula_to_fit, 'StartPoint', sp, 'Lower', lb, 'Upper', ub);
    end
end

function [y,h0] = plotData(x,m,se,f,chan_plot_idx,colors,fits_only,c)  
    if ~exist('c', 'var'), c = zeros(1,length(chan_plot_idx)); end % constants to add to each channel
    ally = nan(length(x), length(chan_plot_idx)); % to use for scaling y-axis
    for ii = 1:length(chan_plot_idx)
        y = m(:,chan_plot_idx(ii))+ c(ii);
        % plot data
        if ~fits_only
            if ~isempty(se)
                y_se = se;
                errorbar(x, y, y-y_se(:,1), y_se(:,2)-y, '.', 'Color', colors(ii,:), 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0);
            else
                plot(x,y, '.','Color', colors(ii,:), 'MarkerSize', 30, 'LineStyle', 'none');
            end
        end
        ally(:,ii) = y;
        % plot fits
        if ~isempty(f)
            f1 = f{ii};
            x1 = 0:max(x)/1000:max(x);
            y1 = f1(x1) + c(ii);
            %h0 = plot(x1, y1, 'Color', colors(ii,:), 'LineWidth', 2); 
%             if ~fits_only
%                 set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%             end
        end
    end
    % format axes
    y = ally(:);
    set(gca, 'XTick', x([1 end]),'XLim', [0 max(x)+0.1*max(x)], 'YLim', [0 max(y)+0.1*max(y)], 'XTickLabelRotation', 0);
end
