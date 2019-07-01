


% Pre-processed file locations for subject 648 and subject 661:
% /Volumes/server/Projects/BAIR/ECoG/648/NY648_data_epoched.mat

% stimulus numbers correspond to the temporal experiment: 25 - 36 (25 - 30, one pulse; 31 - 36 two pulses)

%% LOAD DATA

%dataLoc = '/Volumes/server/Projects/BAIR/ECoG/648/';
dataLoc = '/Volumes/server/Projects/Temporal_integration/DN2_2018_code/data/';
dataNm  = 'NY648_data_epoched_bbmethod5.mat';
a       = load(fullfile(dataLoc, dataNm));

%% PRE-DEFINED PARAMETERS

data = [];

stimNm = [1 : 5, 25 : 36];
stimLb = a.trials.stimuli.soc;

eltomatch = {'MO_01', 'MO_02', 'MO_03', 'MO_04'}; % electrodes with coverage in v1-v3
el        = ecog_matchchannels(eltomatch, a.trials);
nElec     = length(el);

% BIN DATA ACCORDING TO TEMPORAL CONDITIONS
for k = 1 : length(stimNm)
    data(:, :, k) = squeeze(mean(a.trials.broadband.soc(el, :, stimLb == stimNm(k)), 3)); % dimensions: channel x time x stimuli
end


% BIN DATA ACCORDING TO TEMPORAL CONDITIONS
for k = 1 : length(stimNm)
    data(:, :, k) = squeeze(median(a.trials.broadband.soc(el, :, stimLb == stimNm(k)), 3)); % dimensions: channel x time x stimuli
end

mdata = squeeze(mean(data));

%% VISUALIZE DATA

% PLOT THE CONTRAST TIME COURSES
figure (1), clf
for k = 1 : nElec
    subplot(3, nElec, k), set(gca, 'colororder', copper(5)), hold on
    plot(squeeze(data(k, :, 1 : 5)), 'linewidth', 2), axis tight, xlim([200, 700])
    
    subplot(3, nElec, k + nElec), set(gca, 'colororder', copper(6)), hold on
    plot(squeeze(data(k, :, 6 : 11)), 'linewidth', 2), axis tight, xlim([200, 700])
    
    subplot(3, nElec, k + 2 * nElec), set(gca, 'colororder', copper(6)), hold on
    plot(squeeze(data(k, :, 12 : 17)), 'linewidth', 2), axis tight, xlim([200, 700])
end

%% THE SECOND SUBJECT


%dataLoc2 = '/Volumes/server/Projects/BAIR/ECoG/661/';

dataNm2  = 'NY661_data_epoched_bbmethod4.mat';
b = load(fullfile(dataLoc, dataNm2));

%% PRE-DEFINED PARAMETERS FOR THE SECOND SUBJECT

data = [];

stimNm = [1 : 5, 25 : 36];
stimLb = b.trials.stimuli.soc;

eltomatch = {'OC_01', 'OC_02', 'OC_05', 'DPMT_01'}; % electrodes with visual coverage
el        = ecog_matchchannels(eltomatch, b.trials);
nElec     = length(el);

% BIN DATA ACCORDING TO TEMPORAL CONDITIONS
for k = 1 : length(stimNm)
    data(:, :, k) = squeeze(mean(b.trials.broadband.soc(el, :, stimLb == stimNm(k)), 3)); % dimensions: channel x time x stimuli
end

for k = 1 : length(stimNm)
    data(:, :, k) = squeeze(median(b.trials.broadband.soc(el, :, stimLb == stimNm(k)), 3)); % dimensions: channel x time x stimuli
end

%% VISUALIZE DATA

% PLOT THE CONTRAST TIME COURSES
figure (3), clf
for k = 1 : nElec
    subplot(3, nElec, k), set(gca, 'colororder', copper(5)), hold on
    plot(squeeze(data(k, :, 1 : 5)), 'linewidth', 2), axis tight, xlim([200, 700])
    
    subplot(3, nElec, k + nElec), set(gca, 'colororder', copper(6)), hold on
    plot(squeeze(data(k, :, 6 : 11)), 'linewidth', 2), axis tight, xlim([200, 800])
    
    subplot(3, nElec, k + 2 * nElec), set(gca, 'colororder', copper(6)), hold on
    plot(squeeze(data(k, :, 12 : 17)), 'linewidth', 2), axis tight, xlim([200, 800])
end


