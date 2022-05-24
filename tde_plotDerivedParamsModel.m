function tde_plotDerivedParamsModel(params,objFunction,channels,t,chan_to_plot,fits_only,savePlot,saveStr)

% For each channel in the data, plot 6 things:
% 1. contrast response function (based on CRF trials)
% 2. contrast peak latency (based on CRF trials)
% 3. contrast transient to sustained ratio (based on CRF trials)
% 4. subadditive temporal responses (based on ONEPULSE trials)
% 5. subadditive temporal responses (based on TWOPULSE trials)
% 6. time to recovery (based on TWOPULSE trials)
% 
% 2020 Iris Groen

if ~exist('chan_to_plot', 'var'), chan_to_plot = []; end
if ~exist('fits_only', 'var') || isempty(fits_only), fits_only = false; end
if ~exist('savePlot', 'var') || isempty(savePlot), savePlot = false; end
if ~exist('saveStr', 'var'), saveStr = datestr(now,30); end

% Generate new data based on params
nStim = 50;
[stim, stim_info] = tde_simulateNewStimuli(t,nStim);

% Predict responses for new stimuli
srate = channels.sampling_frequency(1);
pred = nan([size(stim) height(channels)]);
for ii = 1:size(params,2)
    prm = params(:,ii);
   	[~, pred(:,:,ii)] = objFunction(prm, [], stim, srate);      
end
data = pred;

% Determine if data was averaged across elecs prior to fit
if isfield(summary(channels), 'number_of_elecs')
    dataWasAveraged = true;
	figureName = sprintf('summary_electrodeaverages_%s', saveStr);
else
    dataWasAveraged = false;
    %[chan_idx, channels] = groupElecsByVisualArea(channels, 'probabilisticresample', {'V1', 'V2', 'V3', 'V3ab', 'LOTO', 'IPS'});   
    [chan_idx, channels] = groupElecsByVisualArea(channels, 'probabilisticresample', {'V1', 'V2', 'V3', 'higher'});   
    %[chan_idx, channels] = groupElecsByVisualArea(channels, 'fixedassignment', {'V1','V2', 'V3'});   
    %[chan_idx, channels] = groupElecsByVisualArea(channels, 'probabilisticresample', {'V1'});   
	figureName = sprintf('summary_individualelecs_%s', saveStr);
end

% Determine if we're plotting all channels or just a subset
if ~isempty(chan_to_plot) 
    if ~iscell(chan_to_plot), chan_to_plot = {chan_to_plot};end
    chan_plot_idx = nan(length(chan_to_plot),1);
    for ii = 1:length(chan_to_plot)
        %chan_plot_idx(ii) = ecog_matchChannels(chan_to_plot{ii}, channels.name);
        chan_plot_idx(ii) = find(matchAreaNameToAtlas(chan_to_plot{ii}, channels.name));
    end
    chan_plot_idx = sort(chan_plot_idx, 'descend');
else
    chan_plot_idx = height(channels):-1:1;
end

% Set plotting specs
%figure('Position', [360    44   879   654]); hold on;
figure('Name', figureName, 'Position', get(0, 'ScreenSize')); hold on;
colors = flipud(copper(length(chan_plot_idx))); 

% CONTRAST CONDITION (3 PLOTS)

stim_inx = find(contains(stim_info.name, 'CRF'));
stim_on = t>0.05 & t<0.55;
x = stim_info.contrast(stim_inx) * 100; 

% 1. contrast response function
subplot(2,3,1); hold on

% compute sum across stim_on window
m = squeeze(sum(data(stim_on, stim_inx, :),1)); se = []; 

% plot linear prediction 
h0 = line([0 x(end)], [0 1], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% plot data
if ~dataWasAveraged
    [m, se] = averageWithinArea(m, chan_idx);
end
[m, se] = normalizeData(m,se,length(x));
%formula_to_fit = 'x ./ (a + x) * (a+1)';
% formula_to_fit = 'x ./ (a + x) * (a+max(x))/max(x)';
% sp = 1;
% lb = 0;
% [f] = fitData(x,m,chan_plot_idx,formula_to_fit, sp, lb);
[vq] = interpolateData(x,m,chan_plot_idx);
plotData(x,m,se,vq,chan_plot_idx,colors,fits_only);  

% format axes
xlabel('contrast level'); ylabel('summed broadband (0-0.5s)'); title('contrast: contrast response function'); 
axis square;
legend(channels.name(chan_plot_idx), 'Location', 'SouthEast');

% 2. peak latency
subplot(2,3,2); hold on

% compute peak across trial window
m = squeeze(data(:,stim_inx,:)); se = [];
[~,I] = max(m,[],1);
m = squeeze(t(I));

% plot data
if ~dataWasAveraged 
    [m, se] = averageWithinArea(m, chan_idx);
end
m = m*1000; % convert to ms
se = se * 1000;
% formula_to_fit = 'a * x ^ -b';
% sp = [1 1];
% [f] = fitData(x,m,chan_plot_idx,formula_to_fit, sp);
% plotData(x,m,se,f,chan_plot_idx,colors,fits_only);  
[vq] = interpolateData(x,m,chan_plot_idx);
plotData(x,m,se,vq,chan_plot_idx,colors,fits_only);  

% format axes
l = get(gca, 'YLim'); ylim([50 l(2)]);
xlabel('contrast level'); ylabel('peak latency (ms)'); title('contrast: latency shift'); 
legend off; axis square;

% 3. ratio of sustained to  transient
subplot(2,3,3); hold on
m = squeeze(data(:,stim_inx,:)); % timecourse
%m = m./max(m); % normalize timecourse to peak
[M] = max(m,[],1); % value at peak
t_off = (t == 0.5); % offset timepoint
O = m(t_off,:,:); % value at offset
R = M./O; % divide value at offset with value of peak

% plot data
if ~dataWasAveraged 
    [m, se] = averageWithinArea(squeeze(R), chan_idx);
end
% formula_to_fit = 'a * x ^ b + c';
% sp = [1 1 1];
% lb = [0 -inf 1];
% ub = [inf inf 1];
% [f] = fitData(x,m,chan_plot_idx,formula_to_fit, sp, lb, ub);
% plotData(x,m,se,f,chan_plot_idx,colors,fits_only);  
[vq] = interpolateData(x,m,chan_plot_idx);
plotData(x,m,se,vq,chan_plot_idx,colors,fits_only);  

% format axes
xlabel('contrast level'); ylabel('ratio offset to peak'); title('contrast: transient to sustained ratio'); 
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
h0 = line([0 x(end)], [0 1], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% plot data
if ~dataWasAveraged
    [m, se] = averageWithinArea(m, chan_idx);
end
[m, se] = normalizeData(m,se, length(x));
formula_to_fit = 'x ./ (a + x) * (a+max(x))/max(x)';
% sp = 1;
% lb = 0;
% [f] = fitData(x,m,chan_plot_idx,formula_to_fit, sp, lb);
% plotData(x,m,se,f,chan_plot_idx,colors,fits_only);  
[vq] = interpolateData(x,m,chan_plot_idx);
plotData(x,m,se,vq,chan_plot_idx,colors,fits_only);  

% format axes
xlabel('stimulus duration (ms)'); ylabel('summed broadband (0-0.5s)'); title('duration: temporal summation'); 
legend off; axis square;


% ISI CONDITION PLOTS (2)

stim_inx = find(contains(stim_info.name, {'ONEPULSE-5', 'TWOPULSE'}));
stim_on = t>0 & t(end); % compute sum over entire epoch, so we don't miss the second pulse; 
x = stim_info.ISI(stim_inx) * 1000; % in ms

% % 1. temporal summation
% subplot(2,3,5); hold on
% 
% % compute sum across stim_on window
% m = squeeze(sum(data(stim_on, stim_inx, :),1)); se = []; 
% 
% % subtract off the 0 ISI condition
% % m = m-m(1,:);
% 
% % plot data
% if ~dataWasAveraged
%     [m, se] = averageWithinArea(m, chan_idx);
% end
% 
% % % subtract off the 0 ISI condition
% % c_min = m(1,:);
% % c_max = m(length(x),:);
% % c = c_min./c_max;
% % m = m-c_min;
% % for cc = 1:2
% %     se(:,:,cc) = se(:,:,cc)-c_min;
% % end
% 
% [m, se] = normalizeData(m,se,length(x));
% 
% % formula_to_fit = 'x ./ (a + x) + b ';
% % sp = [1 1];
% % lb = [0 -10];
% % ub = [nan 0];
% % [f] = fitData(x,m,chan_plot_idx,formula_to_fit, sp, lb, ub);
% 
% % formula_to_fit = 'a * x ^ b + c';
% % sp = [1 1 1];
% % [f] = fitData(x,m,chan_plot_idx,formula_to_fit, sp);
% 
% %formula_to_fit = 'x .^ b ./ (max(x) .^ b)';
% formula_to_fit = '(c + (x .^ b)) ./ (c + (max(x) .^ b))';
% sp = [1 0];
% lb = [0 -inf];
% ub = [inf inf];
% [f] = fitData(x,m,chan_plot_idx,formula_to_fit, sp, lb, ub);
% 
% [y] = plotData(x,m,se,f,chan_plot_idx,colors,fits_only);  
% 
% % plot linear prediction 
% %h0 = line([x(1) x(end)], [0 0], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
% set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% 
% % format axes
% xlabel('ISI (ms)'); ylabel('summed broadband (0-1.2s)'); title('ISI: temporal summation'); 
% legend off; axis square;

% 2. response recovery with ISI
subplot(2,3,5); hold on

[m] = tde_computeISIrecovery2(data,t,stim_info);
if ~dataWasAveraged
    [m, se] = averageWithinArea(m, chan_idx);
end
%formula_to_fit = 'a * x ^ b + c';
%sp = [1 1 1];
%[f] = fitData(x,m,chan_plot_idx,formula_to_fit, sp);
%[y] = plotData(x,m,se,f,chan_plot_idx,colors,fits_only);  
[vq] = interpolateData(x,m,chan_plot_idx);
plotData(x,m,se,vq,chan_plot_idx,colors,fits_only);  

% plot linear prediction 
h0 = line([x(1) x(end)], [100 100], 'LineStyle', '--', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% format axes
xlabel('ISI (ms)'); ylabel('% of first pulse'); title('ISI: recovery of second pulse'); 
ylim([0 120]);
axis square; legend off;
%legend(channels.name(chan_plot_idx));

% save Plot?
if savePlot
    saveDir = fullfile(analysisRootPath, 'figures', 'summary');
    if ~exist(saveDir, 'dir'), mkdir(saveDir);end
    fprintf('[%s] Saving figures to %s \n',mfilename, saveDir);
    saveas(gcf, fullfile(saveDir, figureName), 'png');
end
end

%%% SUBFUNCTIONS %%%

function [m_norm, se_norm] = normalizeData(m,se,cond)
    %m_norm  = m./max(m);
    m_norm = m./m(cond,:);
    if ~isempty(se)
        se_norm = nan(size(se));
        for cc = 1:2
            se_norm(:,:,cc) = se(:,:,cc)./m(cond,:);
        end
    else
        se_norm = [];
    end
end

function [f] = fitData(x,m,chan_plot_idx,formula_to_fit, sp, lb, ub)
    if ~exist('sp', 'var'), sp = []; end
    if ~exist('lb', 'var'), lb = []; end
    if ~exist('ub', 'var'), ub = []; end
    for ii = 1:length(chan_plot_idx)
        y = m(:,chan_plot_idx(ii));
        f{ii} = fit(x, y, formula_to_fit, 'StartPoint', sp, 'Lower', lb, 'Upper', ub);
    end
end

function [vq] = interpolateData(x,m,chan_plot_idx, method)
	if ~exist('method', 'var'), method = 'linear'; end
    xq = 0:max(x)/1000:max(x);
    for ii = 1:length(chan_plot_idx)
        y = m(:,chan_plot_idx(ii));
        vq{ii} = interp1(x,y,xq, method);
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
                y_se = se(:,chan_plot_idx(ii),:);
                errorbar(x, y, y-y_se(:,:,1), y_se(:,:,2)-y, '.', 'Color', colors(ii,:), 'MarkerSize', 10, 'LineWidth', 1, 'LineStyle', 'none', 'CapSize', 0);
            else
                plot(x,y, '.','Color', colors(ii,:), 'MarkerSize', 10, 'LineStyle', 'none');
            end
        end
        ally(:,ii) = y;
        % plot fits
        if ~isempty(f)
            y1 = f{ii};
            x1 = 0:max(x)/1000:max(x);
            switch class(y1)
                case 'cfit'
                    y1 = y1(x1) + c(ii);
            end
            h0 = plot(x1, y1, 'Color', colors(ii,:), 'LineWidth', 2); 
            if ~fits_only
                set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
        end
    end
    % format axes
    y = ally(:);
    set(gca, 'XTick', x,'XLim', [0 max(x)+0.1*max(x)], 'YLim', [0 max(y)+0.1*max(y)], 'XTickLabelRotation', 90, 'FontSize', 8);
end
