%
%% ADD PATH

addpath(genpath('/Volumes/server/Projects/BAIR/ECoG/code'))

%% LOAD DATA

dataLoc = '/Volumes/server/Projects/BAIR/Analyses/visual/sub-beilen/';
dataNm  = 'sub-beilen_ses-day03_epoched_JYZ.mat';

a = load(fullfile(dataLoc, dataNm));

%% PICK RELEVANT ELECTRODES
n 
el = [];
eltomatch = {'OT15', 'OT16'}; labels = a.trials.label;

for k = 1 : length(eltomatch), el(k) = strmatch(eltomatch{k}, labels); end

%% PRE-DEFINED VARIABLES

% INDEX FOR THE CONTRAST STIMULI IN trials.events.bairspatialpattern:
ctrIdx = [1 : 3; 4: 6; 7 : 9; 10 : 12; 13 : 15]';

% INDEX FOR THE CONTRAST STIMULI IN trials.events.bairtemporalpattern:
for k = 1 : 12, durIdx(k, :) = 1 + (k - 1) * 3 : k * 3; end
durIdx = durIdx';

t = a.trials.time;

%% CHANGE SIGNAL UNIT TO PERCENTAGE CHANGE

% DATA FORMAT: electrodes x time points x trials x runs

base_t  = t >= -0.2 & t < 0;

data{1} = a.trials.broadband.bairspatialpattern(el, :, :, :);
data{2} = a.trials.broadband.bairtemporalpattern(el, :, :, :);

% compute baseline
for iCond = 1 : 2 % spatial  or tempoarl data
    for iElec = 1 : 2
        for irun = 1 : 2
            for iStim = 1 : 36
                tmp = squeeze(data{iCond}(iElec, :, iStim, irun));
                baseline = mean(tmp(base_t));
                data{iCond}(iElec, :, iStim, irun) = data{iCond}(iElec, :, iStim, irun)./baseline - 1;
            end
        end
    end
end

for iCond = 1 : 2
    data{iCond} = mean(data{iCond}, 4);
end
%% EXTRACT DATA

mctrDt = []; mdurDt = [];

for k = 1 : size(ctrIdx, 2)
    mctrDt(:, k, :) = mean(data{1}(:, :, ctrIdx(:, k)), 3);
end

for k = 1 : size(durIdx, 2)
    mdurDt(:, k, :) = mean(data{2}(:, :, durIdx(:, k)), 3);
end


% VISUALIZE THE AVERAGED DATA (1) -----------------------------------------

figure (1), clf,
for k = 1 : 2
    % plot the contrast time courses
    subplot(3, 2, k), set(gca, 'colororder', copper(5)), hold on
    plot(t, squeeze(mctrDt(k, :, :)), 'linewidth', 2), xlim([-0.2, 1]),  ylim([-0.5, 3.2])
    set(gca, 'fontsize', 14), title(sprintf('Electrode %d', k)), text(0.6, 2, 'contrast', 'fontsize', 15)
    % plot the changing duration conditions
    subplot(3, 2, k + 2), set(gca, 'colororder', copper(6)), hold on
    plot(t, squeeze(mdurDt(k, 1 : 6, : )), 'linewidth', 2),  xlim([-0.2, 1]),  ylim([-0.5, 3.2])
    set(gca, 'fontsize', 14), text(0.6, 2, 'Increasing Duration', 'fontsize', 15)
    % plot the changing ISI conditions
    subplot(3, 2, k + 4), set(gca, 'colororder', copper(6)), hold on
    plot(t, squeeze(mdurDt(k, 7 : 12, : )), 'linewidth', 2), xlim([-0.2, 1]), ylim([-0.5, 3.2])
    set(gca, 'fontsize', 14), xlabel('time (s)'), text(0.6, 2, 'Increasing ISI', 'fontsize', 15)
end


