% tde_mkFigure2

% Load data and fits

% electrode-averaged data and DN model fits
d1 = load('/Users/iiagroen/surfdrive/BAIR/Papers/TemporalDynamicsECoG/results/LINEAR_RECTF_EXP_NORM_DELAY_xvalmode0_electrodeaverages_20201030T113053.mat');

% individual electrodes and DN model fits
d2 = load('/Users/iiagroen/surfdrive/BAIR/Papers/TemporalDynamicsECoG/results/LINEAR_xvalmode0_individualelecs_20201029T024020.mat');

% stimulus info
[stim_ts, stim_info] = tde_generateStimulusTimecourses(d1.options.stimnames,t);

% Panel A: example of compressive temporal summation 

conditionsOfInterest = {'ONEPULSE-1', 'ONEPULSE-2'};
timepointsOfInterest = [-0.05 0.35];

stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);


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

figure(1); clf
%set(gcf, 'position',  get(0, 'screensize'));

pos = [0.1 0.6 0.3 0.3];
subplot('position', pos); hold on

plot(s(:), 'color', [0.5 0.5 0.5], 'lineWidth', 1);
plot((d(:)./maxresp), 'k', 'lineWidth', 2);
plot((d_copy(:)./maxresp), 'k:', 'lineWidth', 2);
set(gca, 'xtick',1:size(d,1):length(find(stim_idx))*size(d,1), 'xticklabel', []);
set(gca, 'xticklabel', stim_info.duration(stim_idx)); xlabel('duration (s)'), box off
%title('sub-additive temporal summation');

% Panel B: data and fits
conditionsOfInterest = {'ONEPULSE'};
timepointsOfInterest = [-0.1 0.8];

stim_idx = contains(stim_info.name, conditionsOfInterest);
t_idx    = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

s = stim_ts(t_idx,stim_idx);
d = d1.data(t_idx,stim_idx,1);
p = d1.pred(t_idx,stim_idx,1);

pos = [0.1 0.1 0.8 0.3];
subplot('position', pos); hold on
plot(s(:), 'color', [0.7 0.7 0.7], 'linewidth', 1);
plot(d(:), 'k-', 'linewidth', 2);
plot(p(:), 'r-', 'linewidth', 2);
set(gca, 'xtick',1:size(d,1):size(d,2)*size(d,1), 'xticklabel', stim_info.duration(stim_idx));
box off,  axis tight


% Panel C: temporal summation across electrodes



stim_inx = find(contains(stim_info.name, 'ONEPULSE'));
stim_on = t>0 & t<1;
x = stim_info.duration(stim_inx) * 1000; % in ms

pos = [0.5 0.6 0.4 0.3];
subplot('position', pos); hold on

% compute sum across stim_on window
m = squeeze(sum(d2.data(stim_on, stim_inx, :),1)); se = []; 

% plot linear prediction 
h0 = line([0 x(end)], [0 1], 'LineStyle', ':', 'LineWidth', 2, 'Color', [0.7 0.7 0.7]);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca, 'xtick', x);
%%
% plot data
if ~dataWasAveraged
    [m, se] = averageWithinArea(m, group_prob);
end
[m, se] = normalizeData(m,se, length(x));
formula_to_fit = 'x ./ (a + x) * (a+max(x))/max(x)';
sp = 1;
lb = 0;
[f] = fitData(x,m,chan_plot_idx,formula_to_fit, sp, lb);
plotData(x,m,se,f,chan_plot_idx,colors,fits_only);  

% format axes
xlabel('stimulus duration (ms)'); ylabel('summed broadband (0-0.5s)'); title('duration: temporal summation'); 
legend off; axis square;


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