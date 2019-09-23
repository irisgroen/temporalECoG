

load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som648/ses-nyuecog01/sub-som648_ses-nyuecog01_epoched_20190625T161823.mat'); % old

load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som648/ses-nyuecog01/sub-som648_ses-nyuecog01_epoched_20190917T145750.mat'); % new up to 200 Hz

%load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som648/ses-nyuecog01/sub-som648_ses-nyuecog01_epoched_20190917T150828.mat'); % new up to 170 Hz

load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som648/ses-nyuecog01/sub-som648_ses-nyuecog01_epoched_20190917T153335.mat'); % old bbmethod 170 Hz

load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som648/ses-nyuecog01/sub-som648_ses-nyuecog01_epoched_20190917T160907.mat'); % new w old bbmethod

load('/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed/sub-som648/ses-nyuecog01/sub-som648_ses-nyuecog01_epoched_20190917T165952.mat'); % new w old bbmethod 200

%%
specs = [];
specs.baselineType = 'selectedtrials';
specs.dataTypes = {'broadband'};

whichElectrodes = {'MO01', 'MO02', 'MO03', 'MO04'};
trialType = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
ecog_plotTimecourses(trials, whichElectrodes, trialType, specs)