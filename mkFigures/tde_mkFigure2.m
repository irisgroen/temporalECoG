% tde_mkFigure2

% Load data and fits

% electrode-averaged data and DN model fits
modelfun = @LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datatype = 'electrodeaverages';
[d1] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% individual electrodes and DN model fits
datatype = 'individualelecs';
[d2] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% Load stimulus info
[stim_ts, stim_info] = tde_generateStimulusTimecourses(d1.options.stimnames,d1.t);
stim_info.duration = stim_info.duration*1000; % convert to ms

% Subplot positions: % [left bottom width height]
posa = [0.1 0.55 0.4 0.3];
posb = [0.1 0.1 0.8 0.25];
posc = [0.55 0.45 0.35 0.45];

figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

%% Panel A: example of compressive temporal summation 

% Select two conditions to plot
conditionsOfInterest = {'ONEPULSE-1', 'ONEPULSE-2'};
timepointsOfInterest = [-0.05 0.3];

stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = d1.t>timepointsOfInterest(1) & d1.t<=timepointsOfInterest(2);


s = stim_ts(t_idx,stim_idx);
d = d1.data(t_idx,stim_idx,1);
sft = length(find(s(:,1)));
d_shift = padarray(d(:,1), [sft, 0], 0, 'pre');
s_shift = padarray(s(:,1), [sft, 0], 0, 'pre');
d_shift = d_shift(1:size(s,1));
s_shift = s_shift(1:size(s,1));
d_sum = sum([d_shift d(:,1)],2);
d_copy = d;
d_copy(:,1) = nan;
d_copy(:,2) = d_sum;
maxresp = max(d(:,1)); % scale stimulus to max of lowest duration

% for plotting purposes, add a little space after the first stimulus
dummy = nan(40,2);
s = cat(1, s, dummy);
d = cat(1, d, dummy);
d_copy = cat(1, d_copy, dummy);

subplot('position', posa); hold on

plot(s(:), 'color', [0.5 0.5 0.5], 'lineWidth', 1);
plot((d(:)./maxresp), 'k', 'lineWidth', 2);
plot((d_copy(:)./maxresp), 'k:', 'lineWidth', 2);
set(gca, 'xtick',1:size(d,1):length(find(stim_idx))*size(d,1), 'ytick', []);
set(gca, 'xticklabel', stim_info.duration(stim_idx)); box off
xlabel('Stimulus duration (ms)'); ylabel('Neural response'); title('Temporal summation is sub-additive', 'fontsize', 20); 
%axis tight
legend({'Stimulus', 'Neural data', 'Linear prediction'}, 'location', 'northwest', 'fontsize', 18);
legend('boxoff')

%% Panel B: data and fits
conditionsOfInterest = {'ONEPULSE'};
timepointsOfInterest = [-0.1 0.8];

stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = d1.t>timepointsOfInterest(1) & d1.t<=timepointsOfInterest(2);

s = stim_ts(t_idx,stim_idx);
d = d1.data(t_idx,stim_idx,1);
p = d1.pred(t_idx,stim_idx,1);
maxresp = max(d(:,1)); % scale stimulus to max of first condition

subplot('position', posb); hold on
hs = plot(s(:)*maxresp, 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(d(:), 'k-', 'linewidth', 2);
hp = plot(p(:), 'r-', 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xtick',1:size(d,1):size(d,2)*size(d,1), 'xticklabel', stim_info.duration(stim_idx));
box off,  axis tight
xlabel('Stimulus duration (ms)'); ylabel('Change in power (x-fold)'); title('Broadband responses to increasing durations', 'fontsize', 20); 
legend({'Neural data', 'DN model prediction'}, 'location', 'northwest', 'fontsize', 18);
legend('boxoff');

%% Panel C: temporal summation across electrodes
stim_inx = find(contains(stim_info.name, 'ONEPULSE'));
stim_on = d2.t>0 & d2.t<1;
x = stim_info.duration(stim_inx); % in ms

subplot('position', posc); hold on

% plot linear prediction 
h0 = line([0 x(end)], [0 1], 'linestyle', ':', 'LineWidth', 2, 'color', [0 0 0]);
%set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca, 'xtick', x);

% Compute sum across stim_on window
m = squeeze(sum(d2.data(stim_on, stim_inx, :),1)); 
[~, channels, group_prob] = groupElecsByVisualArea(d2.channels, 'probabilisticresample', {'V1'});   

% Plot average parameter values within groups
[m, se] = averageWithinArea(m, group_prob, [], 10000);
se(:,1) = se(:,1)./m(end);
se(:,2) = se(:,2)./m(end);
m = m./m(end);
formula_to_fit = 'x ./ (a + x) * (a+max(x))/max(x)';
sp = 1;
lb = 0;
[f] = fitData(x,m,1,formula_to_fit, sp, lb);
plotData(x,m,se,f,1,[0 0 0],0);  

% Generate new stimuli based on params
nStim = 50;
[stim, stim_info] = tde_simulateNewStimuli(d2.t,nStim);
stim_inx = find(contains(stim_info.name, 'ONEPULSE'));
stim_on = d2.t>0 & d2.t<1;
x = stim_info.duration(stim_inx); % in ms

% Predict model responses for new stimuli
srate = d2.channels.sampling_frequency(1);
pred = nan([size(stim(:,stim_inx)) height(d2.channels)]);
for ii = 1:size(d2.params,2)
    prm = d2.params(:,ii);
   	[~, pred(:,:,ii)] = d2.objFunction(prm, [], stim(:,stim_inx), srate);      
end

% Plot average model predictions across electrodes
m = squeeze(sum(pred(stim_on,:,:),1)); 
[m, se] = averageWithinArea(m, group_prob, [], 10000);
se(:,1) = se(:,1)./m(end);
se(:,2) = se(:,2)./m(end);
m = m./m(end);
x = stim_info.duration(stim_inx)*1000;
plot(x,m, 'r', 'linewidth', 2);
ch = ciplot(se(:,1), se(:,2), x, 'r', 0.25);
set(get(get(ch,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% Format axes
xlabel('Stimulus duration (ms)'); ylabel('Summed broadband timecourse (0-1s)'); title('Temporal summation', 'fontsize', 20); 
legend({'Linear prediction', 'Neural data', 'DN model prediction'}, 'location', 'southeast', 'fontsize', 18);

legend('boxoff')
axis square
axis tight

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
