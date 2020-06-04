function tde_plotDerivedParamsData(data,channels,t,stim_info, chan_to_plot)

% For each channel in the data, compute 4 things:
% - contrast response function (based on CRF trials)
% - subadditive temporal responses (based on ONEPULSE trials)
% - subadditive temporal responses (based on TWOPULSE trials)
% - time to recovery (based on TWOPULSE trials)

if ~exist('chan_to_plot', 'var'), chan_to_plot = []; end

[nSamp, nStim, nChan] = size(data);

% Determine if data was averaged across elecs prior to fit
if isfield(summary(channels), 'number_of_elecs')
    dataWasAveraged = true;
else
    dataWasAveraged = false;
    [chan_idx, channels] = groupElecsByVisualArea(channels, 'probabilisticresample');   
    nChan = height(channels);
end

% Determine if we're plotting all channels or just a subset
if ~isempty(chan_to_plot) 
    if ~iscell(chan_to_plot), chan_to_plot = {chan_to_plot};end
    chan_plot_idx = nan(length(chan_to_plot));
    for ii = 1:length(chan_to_plot)
        chan_plot_idx(ii) = ecog_matchChannels(chan_to_plot{ii}, channels.name);
        nChan = length(chan_plot_idx);
    end 
else
    chan_plot_idx = nChan:-1:1;
end

% Set plotting specs
figure('Position', [360    44   879   654]); hold on;
colors = flipud(copper(length(chan_plot_idx))); 

% CONTRAST CONDITION PLOTS (3)

stim_inx = find(contains(stim_info.name, 'CRF'));
stim_on = t>0 & t<0.55;
x = stim_info.contrast(stim_inx) * 100; 

% 1. contrast response function
subplot(2,3,1); hold on

% compute sum across stim_on window
m = squeeze(sum(data(stim_on, stim_inx, :),1)); se = []; 

% plot linear prediction 
h0 = line([x(1) x(end)], [0 1], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% plot data
if ~dataWasAveraged
    [m, se] = averageWithinArea(m, chan_idx);
end
[m, se] = normalizeData(m,se);
formula_to_fit = 'x ./ (a + x)';
startpoint = [1];
[y] = plotData(x,m,se,chan_plot_idx,colors,formula_to_fit, startpoint);

% format axes
xlabel('contrast level'); ylabel('summed broadband (0-0.5s)'); title('contrast response functions'); 
legend off; axis square;

% 2. peak latency
subplot(2,3,2); hold on

% compute peak across trial window
m = squeeze(data(:,stim_inx,:)); se = [];
[~,I] = max(m,[],1);

% plot data
if ~dataWasAveraged 
    [t_peak, se] = averageWithinArea(squeeze(t(I)), chan_idx);
end
t_peak = t_peak*1000; % convert to ms
se = se * 1000;
formula_to_fit = 'a * x ^ -b';
startpoint = [1 1];
[y] = plotData(x,t_peak,se,chan_plot_idx,colors,formula_to_fit, startpoint);

% format axes
l = get(gca, 'YLim'); ylim([50 l(2)]);
xlabel('contrast level'); ylabel('peak latency (ms)'); title('latency shift'); 
legend off; axis square;

% 3. ratio of offset to peak
subplot(2,3,3); hold on
m = squeeze(data(:,stim_inx,:));
m = m./max(m);
[M] = max(m,[],1);
t_off = t == 0.5;
O = m(t_off,:,:);
R = M./O;

% plot data
if ~dataWasAveraged 
    [m, se] = averageWithinArea(squeeze(R), chan_idx);
end
formula_to_fit = 'a * x ^ b + c';
startpoint = [1 1 1];
[y] = plotData(x,m,se,chan_plot_idx,colors,formula_to_fit, startpoint);

% format axes
xlabel('contrast level'); ylabel('ratio offset to peak'); title('transient to sustained ratio'); 
legend off; axis square;

% DURATION CONDITION PLOTS (1)

stim_inx = find(contains(stim_info.name, 'ONEPULSE'));
stim_on = t>0 & t<0.55;
x = stim_info.duration(stim_inx) * 1000; % in ms

% 1. temporal summation
subplot(2,3,4); hold on

% compute sum across stim_on window
m = squeeze(sum(data(stim_on, stim_inx, :),1)); se = []; 

% plot linear prediction 
h0 = line([x(1) x(end)], [0 1], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% plot data
if ~dataWasAveraged
    [m, se] = averageWithinArea(m, chan_idx);
end
[m, se] = normalizeData(m,se);
formula_to_fit = 'x ./ (a + x)';
startpoint = 0;
[y] = plotData(x,m,se,chan_plot_idx,colors,formula_to_fit, startpoint);

% format axes
xlabel('stimulus duration (ms)'); ylabel('summed broadband (0-0.5s)'); title('temporal summation with duration'); 
legend off; axis square;


% ISI CONDITION PLOTS (2)

% 1. temporal summation
stim_inx = find(contains(stim_info.name, {'ONEPULSE-5', 'TWOPULSE'}));
stim_on = t>0 & t(end); % compute sum over entire epoch, we don't miss the second pulse; 
x = stim_info.ISI(stim_inx) * 1000; % in ms

% 1. temporal summation
subplot(2,3,5); hold on

% compute sum across stim_on window
m = squeeze(sum(data(stim_on, stim_inx, :),1)); se = []; 

% plot data
if ~dataWasAveraged
    [m, se] = averageWithinArea(m, chan_idx);
end
[m, se] = normalizeData(m,se);
formula_to_fit = 'b + x ./ (a + x)';
startpoint = [1 1];
[y] = plotData(x,m,se,chan_plot_idx,colors,formula_to_fit, startpoint);

% plot linear prediction 
h0 = line([x(1) x(end)], [y(1) y(1)], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% format axes
xlabel('ISI (ms)'); ylabel('summed broadband (0-1.2s)'); title('temporal summation with ISI'); 
legend off; axis square;

% 2. response recovery with ISI
subplot(2,3,6); hold on

[m] = tde_computeISIrecovery(data,t,stim_info);
if ~dataWasAveraged
    [m, se] = averageWithinArea(m, chan_idx, @mean);
end
formula_to_fit = 'a * x ^ b + c';
startpoint = [1 1 1];
[y] = plotData(x,m,se,chan_plot_idx,colors,formula_to_fit,startpoint);

% plot linear prediction 
h0 = line([x(1) x(end)], [100 100], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% format axes
xlabel('ISI (ms)'); ylabel('% of first pulse'); title('recovery with ISI'); 
ylim([0 100]);
legend(channels.name(chan_plot_idx));

end

%%% SUBFUNCTIONS %%%

function [m_norm, se_norm] = normalizeData(m,se)
    m_norm  = m./max(m);
    if ~isempty(se)
        se_norm = nan(size(se));
        for cc = 1:2
            se_norm(:,:,cc) = se(:,:,cc)./max(m);
        end
    else
        se_norm = [];
    end
end


function [y,h0] = plotData(x,m,se,chan_plot_idx,colors,formula_to_fit, startpoint)
if ~exist('startpoint', 'var'), startpoint = []; end
ally = nan(length(x), length(chan_plot_idx));
for ii = 1:length(chan_plot_idx)
    y = m(:,chan_plot_idx(ii));
    if ~isempty(se)
        y_se = se(:,chan_plot_idx(ii),:);
        errorbar(x, y, y-y_se(:,:,1), y_se(:,:,2)-y, '.', 'Color', colors(ii,:), 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0);
    else
        plot(x,y, '.','Color', colors(ii,:), 'MarkerSize', 30, 'LineStyle', 'none');
    end
    ally(:,ii) = y;
    if ~isempty(formula_to_fit)
        f = fit(x, y, formula_to_fit, 'StartPoint', startpoint);
        if length(chan_plot_idx) > 1
            h0 = plot(x,f(x), 'Color', colors(ii,:), 'LineWidth', 2); 
        else
            h0 = plot(f);
        end
        set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
end
y = ally;
% format axes
set(gca, 'XTick', x,'XLim', [0 max(x)+0.1*max(x)], 'YLim', [0 max(y(:))+0.1*max(y(:))], 'XTickLabelRotation', 45);

end
