% fit the DN model to the new NYU data

%% LOAD DATA
subj = 648;

addpath(genpath('/Volumes/server/Projects/Temporal_integration/DN_2018_code_data/code'))

dataLoc = fullfile(dn_ctrst_RootPath, 'data');

dataNm  = sprintf('NY%d_data_epoched.mat', subj);
a       = load(fullfile(dataLoc, dataNm));

%% PRE-DEFINED VARIABLES

data = [];

stimNm = [1 : 5, 25 : 36];
stimLb = a.trials.stimuli.soc;

switch subj
    case 661
        eltomatch = {'OC_01', 'OC_02', 'OC_05', 'DPMT_01'};
    case 648
        eltomatch = {'MO_01', 'MO_02', 'MO_03', 'MO_04'}; % electrodes with coverage in v1-v3
end

el        = ecog_matchchannels(eltomatch, a.trials);
nElec     = length(el);

contrasts = [0.0405, 0.0902, 0.2105, 0.3203, 1.0000];

durs =[0.016667,0.033333, 0.066667, 0.13333, 0.26667, 0.53333];

t = a.trials.time;

%% CHANGE BROADBAND DATA BASELINE

bb = a.trials.broadband.soc; % number of electrodes, time course, stimulus index

% COMPUTE BASELINE --------------------------------------------------------
% use the firsrt 200ms before stimulus onset as the baseline
base_range = t >= -0.2 & t < 0;

%base_idx = 1 : 200;
m_base   = squeeze(mean(mean(bb(:, base_range, :), 2), 3)); % one mean baseline number per electrode

% EXTRACT THE RELEVANT STIMULUS
dataAllTrials = [];
data = [];

for k = 1 : length(stimNm)
    rel_bb = bb(el, :, :)./ m_base(el);
    rel_bb = rel_bb - 1;
    data(:, :, k) = squeeze(mean(rel_bb(:, :, stimLb == stimNm(k)), 3)); % dimensions: channel x time x stimuli
    dataAllTrials(:, :, k,:) = squeeze(rel_bb(:, :, stimLb == stimNm(k))); % dimensions: channel x time x stimuli
end

% NORMALIZE DATA TO ITS MAXIMUM HEIGHT FOR EACH ELECTRODE
% for k = 1 : length(eltomatch)
%     tmp = data(k, :, :);
%     data(k, :, :) = tmp ./max(tmp(:));
% end 

%% PLOT THE CONTRAST TIME COURSES

figure (1+subj), clf;
for k = 1 : nElec
    subplot(3, nElec, k), set(gca, 'colororder', copper(5)), hold on
    plot(t, squeeze(data(k, :, 1 : 5)), 'linewidth', 2), axis tight,  xlim([-0.2, 1])
    
    subplot(3, nElec, k + nElec), set(gca, 'colororder', copper(6)), hold on
    plot(t, squeeze(data(k, :, 6 : 11)), 'linewidth', 2), axis tight, xlim([-0.2, 1])
    
    subplot(3, nElec, k + 2 * nElec), set(gca, 'colororder', copper(6)), hold on
    plot(t, squeeze(data(k, :, 12 : 17)), 'linewidth', 2), axis tight, xlim([-0.2, 1])
end


%%

for elec = 1:2
    figure, set(gcf, 'Name', sprintf('Subject %d, Electrode %d', subj, elec)); 
    
    for stim = 1:17 %: nElec
        if stim < 6, plotIdx = stim; else, plotIdx = stim+1; end
        subplot(3, 6, plotIdx), hold on
        plot(t, squeeze(dataAllTrials(elec, :, stim,:)), 'Color', .7*[1 1 1], 'linewidth', 1),
        plot(t, mean(squeeze(dataAllTrials(elec, :, stim,:)),2), 'linewidth', 2)
        plot(t, median(squeeze(dataAllTrials(elec, :, stim,:)),2), 'linewidth', 2)
        % plot(t, geomean(squeeze(dataAllTrials(elec, :, stim,:)),2), 'linewidth', 2)
        axis tight,  axis([-0.2, 1, [0 80]])
    end
    
end
%% AVERAGE BETWEEN THE FIRST TWO ELECTRODES

mdata = squeeze(mean(data(1 : 2, :, :)));

figure (2), clf, set(gca, 'colororder', copper(5)), hold on
plot(t, mdata(:, 1 : 5), 'linewidth', 2)

%% MAKE STIMULUS

nStim = size(data, 3);
stim  = zeros(nStim, length(t));

% CONTRAST STIMULI --------------------------------------------------------
stim(1 : 5, t > 0 & t<=0.5) = 1;
for k = 1 : 5, stim(k, :) = stim(k, :) .* contrasts(k); end

% INCREASING DURATIONS ----------------------------------------------------
for k = 1 : 6, stim(k + 5, (t>0) & (t <= durs(k))) = 1; end

% INCREASING ISI ----------------------------------------------------------
stim(12 : nStim, t > 0 & t <= durs(4)) = 1;
for k = 1 : 6,
    t_start = durs(4) + durs(k);
    t_end   = durs(4) * 2 + durs(k);
    stim(11 + k, t > t_start & t <= t_end) = 1;
end

% the first 5 stimuli are 500-ms sustained stimulus with 

figure (3), clf, imagesc(stim), colormap

%% FIT THE DN MODEL

t_range = t > -0.2 & t <= 1;

t_cut    = t(t_range);
stim_cut = stim(:, t_range);
mdt_cut  = mdata(t_range, :)';

seed = [0.03, 0.07, 1.5, 0.15, 0.06, 1];
seed = [0.1, 0.1, 3, 0.1, 0.06, 1];
lb   = [0, 0, 0, 0, 0, 0];
ub   = [1, 1, 10, 1, 1, 1];

prm = [];
prm = fminsearchbnd(@(x) dn2_fineFitCtrstDur(x, mdt_cut, t_cut, stim_cut), seed, lb, ub);

%% GENERATE MODEL PREDICTIONS

prm_tofit = [prm(1), 0, prm(2 : end)];

pred = dn_DNmodel(prm_tofit, stim_cut, t_cut);
pred = pred./max(pred(:));

%% VISUALIZE

figure (4), clf
for k = 1 : nStim
   subplot(17, 1, k)
   plot(t_cut, stim_cut(k, :), 'k'), hold on
   plot(t_cut, mdt_cut(k, :), 'r-', 'linewidth', 3), 
   plot(t_cut, pred(k, :), 'k-', 'linewidth', 2), 
   axis tight, box off, axis off, ylim([0, 1])
end

%% ALTERNATIVE VISUALIZATION

idx = {1 : 5, 6 : 11, 12 : 17};

figure (5),clf

for k = 1 : 3
   subplot(1, 3, k), %set(gca, 'colororder', copper(length(idx{k}))), hold on
   plot(t_cut, mdt_cut(idx{k}, :)', 'r-',  'linewidth', 2), hold on
   
   %subplot(2, 3, k + 3), set(gca, 'colororder', copper(length(idx{k}))), hold on
   plot(t_cut, pred(idx{k}, :)','k-', 'linewidth', 2), 
   axis tight, box off
end



    
    
