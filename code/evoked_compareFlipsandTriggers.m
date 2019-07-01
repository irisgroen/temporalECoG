%% load
%a{1} = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som648/ses-nyuecog01/sub-som648_ses-nyuecog01_epoched.mat');
%a{2} = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som648/ses-nyuecog01/sub-som648_ses-nyuecog01_epoched_segmentedOnTriggers.mat');

%a{1} = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som692/ses-nyuecog01/sub-som692_ses-nyuecog01_epoched.mat');
%a{2} = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som692/ses-nyuecog01/sub-som692_ses-nyuecog01_epoched_segmentedOnTriggers.mat');

%a{1} = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som708/ses-nyuecog01/sub-som708_ses-nyuecog01_epoched.mat');
%a{2} = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som708/ses-nyuecog01/sub-som708_ses-nyuecog01_epoched_segmentedOnTriggers.mat');

a{1} = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som723/ses-nyuecog01/sub-som723_ses-nyuecog01_epoched.mat');
a{2} = load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som723/ses-nyuecog01/sub-som723_ses-nyuecog01_epoched_segmentedOnTriggers.mat');

%% baselinecorrect

t = a{1}.trials.time;

% extract evoked activity
ev = {};
for k = 1 : length(a)
    ev{k} = a{k}.trials.evoked;
end

% baseline correction

% use the firsrt 200ms before stimulus onset as the baseline
base_range = (t >= -0.2 & t < 0);

for k = 1 : length(a)
    %m_base{k} = squeeze(median(mean(ev{k}(:, base_range, :), 2), 3));
    %ev{k}     = ev{k}./m_base{k};
    baseline = mean(ev{k}(:,base_range,:),2);
    ev{k} = ev{k} - baseline;
end

%% plot

%stimNm = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5',...
%^'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
%'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};
differences = [a{1}.trials.events.event_sample - a{2}.trials.events.event_sample];
fsample = a{1}.trials.fsample;
figure;hist(differences)
title('flip - triggers'); xlabel('samples');
figure;hist(differences * (1/fsample))
title('flip - triggers'); xlabel('time (s)'); xlim([-0.2 0.2])


stimNm = {'ONEPULSE-1'};
taskNm = 'prf';
%chan = 'MO01'; % 648
%chan = 'BO01'; % 692
%chan = 'O02'; % 708
chan = 'RPT05'; % 723

versiontitles = {'flips', 'triggers'};

for k = 1:2
    test = squeeze(ev{k}(ecog_matchChannels(chan, a{k}.trials),:,find(contains(a{k}.trials.events.trial_name, stimNm))));
    %test = squeeze(ev{k}(ecog_matchChannels(chan, a{k}.trials),:,find(contains(a{k}.trials.events.task_name, taskNm))));
    figure;plot(t, test); title(versiontitles{k});
    hold on; plot(t,mean(test,2), 'k', 'LineWidth', 2);
    %set(gca, 'xlim', [0 0.2]);
end

