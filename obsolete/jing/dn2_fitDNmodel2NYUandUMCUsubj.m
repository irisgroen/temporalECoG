% dn2_fitDNmodel2BothNYUsubjects

saveFigure = 1;
whichElecs = 'V1'; % 'V1' or 'V2';

%% SET PATHS

% Check whether we have the ECoG_utils repository on the path
if ~exist('ecog_matchChannels', 'file')
    tbUse ECoG_utils;
end

%addpath(genpath('/Volumes/server/Projects/Temporal_integration/DN_2018_code_data/code'))
addpath(genpath('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal'));
addpath(genpath('/Volumes/server/Projects/Temporal_integration/DN_2018_code_data/code'))

%% LOAD DATA

% NYU 648 661
dataLoc = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed';

sub_label = [648, 661];
ses_label = 'nyuecog01';

for k = 1 : length(sub_label)
    dataNm  = fullfile(dataLoc, sprintf('sub-som%d', sub_label(k)), sprintf('ses-%s', ses_label),  sprintf('sub-som%d_ses-%s_epoched.mat', sub_label(k), ses_label));
    a{k}    = load(dataNm);
end

% UMCU CHAAM
dataLoc = '/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess';

sub_label = 'chaam';
ses_label = 'UMCUECOGday03';

dataName = fullfile(dataLoc, sprintf('umcu%s_preproc_selectelecs',sub_label));
a{3}    = load(dataName);

%% PRE-DEFINED VARIABLES

data = [];

%stimNm    = [1 : 5, 25 : 36];

stimNm = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5',...
    'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
    'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};

contrasts = [0.0405, 0.0902, 0.2105, 0.3203, 1.0000];
durs      = [0.016667, 0.033333, 0.066667, 0.13333, 0.26667, 0.53333];

eltomatch{1} = {'MO01', 'MO02'};
eltomatch{2} = {'OC01'};
eltomatch{3} = a{3}.trials.channels.name(:);

t      = a{1}.trials.time;

for k = 1 : 3
    stimLb{k} = a{k}.trials.events.trial_name;    
    el{k}     = ecog_matchChannels(eltomatch{k}, a{k}.trials);
    nElec(k)  = length(el{k});
end

%% CHANGE BROADBAND DATA BASELINE
bb = {};
bb{1} = a{1}.trials.broadband(el{1}, :, :);
bb{2} = a{2}.trials.broadband(el{2}, :, :);
bb{3} = a{3}.trials.broadband(el{3}, :, :);

% DEFINE BASELINE RANGE ---------------------------------------------------
% use the firsrt 200ms before stimulus onset as the baseline
base_range = (t >= -0.2 & t < 0);

for k = 1 : 3
    m_base{k} = squeeze(median(mean(bb{k}(:, base_range, :), 2), 3));
    bb{k}     = bb{k}./m_base{k}-1;
end

% HACK downsample (??) UMCU data
nSamp = size(bb{3},2);
bb{3} = bb{3}(:,1:4:nSamp,:);

% SPLIT DATA INTO DIFFERENT STIMULUS TYPES --------------------------------
data = [];
% for the first subject
for k = 1 : 2
    for k1 = 1 : length(stimNm)
        data(k, k1, :) = squeeze(median(bb{1}(k, :, contains(stimLb{1},stimNm{k1})), 3));
    end
end

% for the second subject
for k1 = 1 : length(stimNm)
   data(3, k1, :) = squeeze(median(bb{2}(:, :, contains(stimLb{2},stimNm{k1})), 3));
end

% for the third subject
for k = 1 : size(el{3},2)
    for k1 = 1 : length(stimNm)
        data(k+3, k1, :) = squeeze(median(bb{3}(k, :, contains(stimLb{3},stimNm{k1})), 3));
    end
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

for k = 1 : size(data,1);
    tmp = data(k, :, :);
    maxRsp(k) = max(tmp(:));
    data(k, :, :) = data(k, :, :)./maxRsp(k);
end

switch whichElecs
    case 'V1'
        data = data(4:9,:,:); % V1 ELECS
    case 'V2'
        data = data([1:3 9:11],:,:); % V2/3 ELECS
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
set(2, 'Position', [100 100 500 1300]);

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

%% plot SUM of prediction

stimData = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/stimuli/sub-umcuchaam_ses-UMCUECOGday03_task-bairtemporalpattern_run-2_acq-clinical_events.mat');
stimDurations = unique(stimData.stimulus.duration);
stimISIs = unique(stimData.stimulus.ISI);

timeInd = t>0 &t<1;
sumData = sum(mdata(:,timeInd),2);
sumPred = sum(pred(:,timeInd),2);

figure (4), clf;

subplot(2,2,1); hold on;
plot(stimDurations, sumPred(6:11), 'Marker', '.', 'MarkerSize', 50,'LineStyle', '-', 'LineWidth', 2, 'Color', 'r');
plot(stimDurations, stimDurations*200,'LineStyle', '-','Color','g', 'LineWidth',2);
xlabel('stimulus duration (s)'); ylabel('sum (0-1s)'); title('model prediction');
axis tight, xlim([0, 0.6]), set(gca, 'fontsize', 16), ylim([0, 120])
legend({'DN model', 'linear prediction'}, 'Location', 'NorthWest');

subplot(2,2,2), hold on;
plot(stimISIs, sumPred([10 12:17]), 'Marker', '.', 'MarkerSize', 50,'LineStyle', '-','LineWidth', 2, 'Color', 'r');
plot(stimISIs, ones(length(stimISIs))*60,'LineStyle', '-','Color','g', 'LineWidth',2);
xlabel('stimulus ISI (s)'); ylabel('sum (0-1s)');title('model prediction');
axis tight, xlim([0, 0.6]), set(gca, 'fontsize', 16), ylim([0, 120])

subplot(2,2,3); hold on;
plot(stimDurations, sumData(6:11), 'Marker', '.', 'MarkerSize', 50,'LineStyle', '-', 'LineWidth', 2, 'Color', 'r');
plot(stimDurations, stimDurations*200,'LineStyle', '-','Color','g', 'LineWidth',2);
xlabel('stimulus duration (s)'); ylabel('sum (0-1s)'); title('ecog data');
axis tight, xlim([0, 0.6]), set(gca, 'fontsize', 16), ylim([0, 120])

subplot(2,2,4), hold on;
plot(stimISIs, sumData([10 12:17]), 'Marker', '.', 'MarkerSize', 50,'LineStyle', '-','LineWidth', 2, 'Color', 'r');
plot(stimISIs, ones(length(stimISIs))*60,'LineStyle', '-','Color','g', 'LineWidth',2);
xlabel('stimulus ISI (s)'); ylabel('sum (0-1s)');title('ecog data');
axis tight, xlim([0, 0.6]), set(gca, 'fontsize', 16), ylim([0, 120])

set(gcf, 'Position', [150 100 1000 700]);

%% CALCULATE R2
rvals = [];
for k = 1:17
    rvals(k) = corr(pred(k,:)',mdata(k,:)');  
end
mean(rvals.^2)

%% PLOT SUM vs FMRI (code copied from /Volumes/server/Projects/Temporal_integration/DN_2018_code_data/code/dn_fitDNECoG2fMRI.m'
load('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess/fMRIdataZhouPCB')
mfmri_2fit(2,:) = mean(mfmri_2fit([2 3],:),1); % average V2 and V3 responses
sfmri_2fit(2,:,:) = mean(sfmri_2fit([2 3],:,:),1); % average V2 and V3 responses

roiNm = {'V1', 'V2', 'V3'};

linrsp = [1, 2, 4, 8, 16, 32, 16, 16, 16, 16, 16, 16, 0];
x1 = [0, 1, 2, 4, 8, 16, 32];

switch whichElecs
    case 'V1'
        whichRoi = 1;
        ecogScaleFac = 0.005;
    case 'V2'
        whichRoi = 2;
        ecogScaleFac = 0.0035;
end
fg1 = figure (5); clf
for iroi = whichRoi%nrois
    
    for jj = 1:2
        switch jj 
            case 1, mriOrder = [13, 1 :  6]; ecogOrder = [6:11]; 
            case 2, mriOrder = [5 7:12]; ecogOrder = [10 12:17]; 
        end
     
        subplot(1,2,jj);
    %subplot(2, 3, iroi + 3*(jj-1))
    %figure;hold on;
    
    plot(x1, mfmri_2fit(iroi, mriOrder), 'ko', 'markersize', 9,  'markerfacecolor', 'k'), hold on
    
    for k = 1 : length(mriOrder)
        low = sfmri_2fit(iroi, 1, mriOrder(k));
        high = sfmri_2fit(iroi, 2, mriOrder(k));
        p = plot(x1(k)*[1 1], [low, high], '-', 'linewidth', 3, 'color', 0.8 * ones(1, 3));
        p.Annotation.LegendInformation.IconDisplayStyle = 'off';  
    end
    %xlim([0.5, length(mriOrder) + 0.5]),
    %ylim([-.05, 0.6]), xlim([-2, 34]), box off, set(gca, 'xtick', [0, 4, 16, 32], 'xticklabel', [0, 67, 267, 533])

    title(roiNm{iroi}), set(gca, 'fontsize', 14)

    plot(x1, linscale(iroi).*linrsp(mriOrder), 'g-', 'linewidth', 3)
    plot(x1, ecogPrd1(iroi, mriOrder), 'r-', 'linewidth', 3)
    if jj == 1
        plot(x1(2:end), sumPred(ecogOrder)*ecogScaleFac, 'm-', 'linewidth', 3)
        legend({'BOLD', 'linear pred', 'ECOG fit PCB', 'ECOG fit HBM'}, 'Location', 'NorthWest'); 
    else
        plot(x1, sumPred(ecogOrder)*ecogScaleFac, 'm-', 'linewidth', 3)
    end
    set(gca, 'xaxislocation', 'origin'), xlabel('time (ms)'), ylabel('% BOLD')
    end
end

set(gcf, 'Position', [150 100 1500 700]);

%% SAVE FIGURES

if saveFigure
    switch whichElecs
        case 'V1'
            fg2Nm = 'visualizeModelFit1_NYU+UMCUsubjects_V13elec';
            fg3Nm = 'visualizeModelFit2_NYU+UMCUsubjects_V1elec';
            fg4Nm = 'sumModelFits_NYU+UMCUsubjects_V1elec';
            fg5Nm = 'comparisonWithfMRI_NYU+UMCUsubjects_V1elec';
        case 'V2'
            fg2Nm = 'visualizeModelFit1_NYU+UMCUsubjects_V2V3elec';
            fg3Nm = 'visualizeModelFit2_NYU+UMCUsubjects_V2V3elec';
            fg4Nm = 'sumModelFits_NYU+UMCUsubjects_V2V3elec';
            fg5Nm = 'comparisonWithfMRI_NYU+UMCUsubjects_V2V3elec';
    end
    saveLoc = fullfile(dn_ctrst_RootPath, 'analysisFigures');
    printnice(2, 0, saveLoc, fg2Nm)
    printnice(3, 0, saveLoc, fg3Nm)
    printnice(4, 0, saveLoc, fg4Nm)
    printnice(5, 0, saveLoc, fg5Nm)
end


