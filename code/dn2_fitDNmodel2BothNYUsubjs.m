% dn2_fitDNmodel2BothNYUsubjects

saveFigure = 1;

%% SET PATHS

%addpath(genpath('/Volumes/server/Projects/Temporal_integration/DN_2018_code_data/code'))
addpath(genpath('/Volumes/server/Projects/BAIR/Conference/VSS 2018 VisuoTemporal'));

% Check whether we have the ECoG_utils repository on the path
if ~exist('ecog_matchChannels.m', 'file')
    tbUse ECoG_utils;
end

%dataLoc = fullfile(dn_ctrst_RootPath, 'data');
dataLoc = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed';

%% LOAD DATA

sub_label = [648, 661];
ses_label = 'nyuecog01';

for k = 1 : length(sub_label)
    %dataNm  = sprintf('NY%d_data_epoched.mat', subj(k));
    dataNm  = fullfile(dataLoc, sprintf('sub-som%d', sub_label(k)), sprintf('ses-%s', ses_label),  sprintf('sub-som%d_ses-%s_epoched.mat', sub_label(k), ses_label));

    %a{k}    = load(fullfile(dataLoc, dataNm));
    a{k}    = load(dataNm);
end

%% PRE-DEFINED VARIABLES

data = [];

%stimNm    = [1 : 5, 25 : 36];

stimNm = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5',...
    'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
    'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};

contrasts = [0.0405, 0.0902, 0.2105, 0.3203, 1.0000];
durs      = [0.016667,0.033333, 0.066667, 0.13333, 0.26667, 0.53333];

eltomatch{1} = {'MO01', 'MO02'};
eltomatch{2} = {'OC01'};

for k = 1 : 2
    stimLb{k} = a{k}.trials.events.trial_name;
    t      = a{k}.trials.time;
    el{k}     = ecog_matchChannels(eltomatch{k}, a{k}.trials);
    nElec(k)  = length(el{k});
end

%% CHANGE BROADBAND DATA BASELINE
bb = {};
bb{1} = a{1}.trials.broadband(el{1}, :, :);
bb{2} = a{2}.trials.broadband(el{2}, :, :);

% DEFINE BASELINE RANGE ---------------------------------------------------
% use the firsrt 200ms before stimulus onset as the baseline
base_range = (t >= -0.2 & t < 0);

for k = 1 : 2
    m_base{k} = squeeze(median(mean(bb{k}(:, base_range, :), 2), 3));
    bb{k}     = bb{k}./m_base{k}-1;
end

% SPLIT DATA INTO DIFFERENT STIMULUS TYPES --------------------------------
data = [];
% for the first subject
for k = 1 : 2
    for k1 = 1 : length(stimNm)
        %data(k, k1, :) = squeeze(median(bb{1}(k, :, stimLb{1} == stimNm(k1)), 3));
        data(k, k1, :) = squeeze(median(bb{1}(k, :, contains(stimLb{1},stimNm{k1})), 3));
    end
end
% for the second subject
for k1 = 1 : length(stimNm)
   data(3, k1, :) = squeeze(median(bb{2}(:, :, contains(stimLb{2},stimNm{k1})), 3));
end

%% VISUALIZE DATA

figure (1), clf
for k = 1 : 17
   subplot(3, 17, k)
   plot(squeeze(data(1, k, :))), axis tight, box off,  ylim([-1, 25]),
   
   subplot(3, 17, k + 17)
   plot(squeeze(data(2, k, :))), axis tight, box off,  ylim([-1, 25]),
   
   subplot(3, 17, k + 34)
   plot(squeeze(data(3, k, :))), axis tight, box off,  ylim([-1, 25]),
end

%% COMPUTE THE MAXIMUM RESPONSE IN EACH SUBJECT'S DATA SET

for k = 1 : 3
    tmp = data(k, :, :);
    maxRsp(k) = max(tmp(:));
    data(k, :, :) = data(k, :, :)./maxRsp(k);
end

% average the response across subjects and electrodes

mdata    = squeeze(mean(data));
maxmdata = max(mdata(:));
mdata    = mdata ./maxmdata;

figure (2), clf
for k = 1 : 17
   subplot(17, 1, k)
   plot(t, mdata(k, :), 'r-', 'linewidth', 2), hold on, axis tight, ylim([-0.1, 1]), box off
end

%% MAKE STIMULUS

nStim = size(data, 2);
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

% VISUALIZE THE STIMULI ---------------------------------------------------
figure (2)
for k = 1 : 17
   subplot(17, 1, k)
   plot(t, stim(k, :), 'k-')
end

%% FIT THE DN model

seed = [0.03, 0.07, 1.5, 0.15, 0.06, 1];
%seed = [0.1, 0.1, 3, 0.1, 0.06, 1];
lb   = [0, 0, 0, 0, 0, 0];
ub   = [1, 1, 10, 1, 1, 1];

prm = [];
prm = fminsearchbnd(@(x) dn2_fineFitCtrstDur(x, mdata, t, stim), seed, lb, ub);

%% GENERATE MODEL PREDICTIONS

prm_tofit = [prm(1), 0, prm(2 : end)];

pred = dn_DNmodel(prm_tofit, stim, t);
pred = pred./max(pred(:));

%% VISUALIZE MODEL FIT

figure (2)
for k = 1 : 17
   subplot(17, 1, k)
   plot(t, pred(k, :), 'k-', 'linewidth', 2)
   set(gca, 'ytick', [0, 1], 'xtick', [0, 0.5, 1]), xlim([-0.2, 1])
   if k == 17, xlabel('time (s)'), end
end

%% ALTERNATIVE WAY TO VISUALIZE DATA

figure (3), clf
subplot(2, 3, 1), set(gca, 'colororder', copper(5)), hold on
plot(t, mdata(1 : 5, :), '-', 'linewidth', 2), 
subplot(2, 3, 4), set(gca, 'colororder', copper(5)), hold on
plot(t, pred(1 : 5, :), '-', 'linewidth', 2), 

subplot(2, 3, 2), set(gca, 'colororder', copper(6)), hold on
plot(t, mdata(6 : 11, :), '-', 'linewidth', 2), 

subplot(2, 3, 5), set(gca, 'colororder', copper(6)), hold on
plot(t, pred(6 : 11, :), '-', 'linewidth', 2), 

subplot(2, 3, 3), set(gca, 'colororder', copper(6)), hold on
plot(t, mdata(12 : 17, :), '-', 'linewidth', 2), 
subplot(2, 3, 6), set(gca, 'colororder', copper(6)), hold on
plot(t, pred(12 : 17, :), '-', 'linewidth', 2),

for k = 1 : 6
    subplot(2, 3, k)
    axis tight, xlim([-0.2, 1]), set(gca, 'ytick', [0, 1], 'xtick', [0, 0.5, 1], 'fontsize', 16), , ylim([0, 1])
    xlabel('time (s)')
end

%% SAVE FIGURES

if saveFigure
    fg2Nm = 'visualizeModelFit1_NYUsubjects3Electrodes';
    fg3Nm = 'visualizeModelFit2_NYUsubjects3Electrodes';
    saveLoc = fullfile(dn_ctrst_RootPath, 'analysisFigures');
    printnice(2, 0, saveLoc, fg2Nm)
    printnice(3, 0, saveLoc, fg3Nm)
end