
% SCRIPT DESCRIPTION %

%% Define paths and dataset specs %%

% Dataset specs
projectName = 'visual';
sub_label   = 'chaam'; 
ses_label   = 'UMCUECOGday03';

% Input paths specs
dataPth     = '/Volumes/server/Projects/BAIR/Data/BIDS/';
dataDir     = fullfile(dataPth, projectName, sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label), 'ieeg');
preprocDir  = fullfile(dataPth, projectName, 'derivatives', 'preprocessed', sprintf('sub-%s', sub_label), sprintf('ses-%s', ses_label));
fsDir       = fullfile('/Volumes/server/Freesurfer_subjects/',sprintf('umcu%s',sub_label));

% Check if we have the ECoG_utils repository on the path
if ~exist('ecog_plotTimecourses.m')
    tbUse ECoG_utils;
end

%% Load preprocessed data

% This file contains stimulus-defined epochs of the rereferenced voltages
% (trials.evoked) and the broadband time series (trials.broadband).
% Prestimulus baseline correction has NOT been performed on these epochs
% (happens inside plotting script - but this can be changed)

dataName = fullfile(preprocDir, sprintf('sub-%s_ses-%s_epoched.mat', sub_label, ses_label));
load(dataName);

%% Match electrodes to visual areas

specs = [];
specs.pID           = sub_label; % patient ID number
specs.patientPool   = 'BAIR';
specs.fsDir         = fsDir;
specs.plotmesh      = 'right';
specs.plotelecs     = 'yes';
trials.viselec      = electrode_to_nearest_node(specs, dataDir);

%% Make a preselection of relevant electrodes

eccThreshold = 12; % degrees
bbThreshold  = 10; % relative change to baseline

% Include channels that have a match with a visual regions in either atlas
ElecsToInclude  = unique([trials.viselec.benson14_varea.elec_labels trials.viselec.wang15_mplbl.elec_labels]); 

% Exclude channels marked as bad or those deemed epileptic
% ---> FROM DORA: Electrodes on seizure are OC12, OC13, OC21, OC22.
ElecsToExclude = [trials.channels.name(contains(trials.channels.status, 'bad'))' {'Oc9', 'Oc12', 'Oc13', 'Oc21', 'Oc22'}]; % Oc9 labeled as bad

% Exclude channels with eccentricity larger than stimulus extent
ElecsToExclude = [ElecsToExclude trials.viselec.benson14_varea.elec_labels(trials.viselec.benson14_varea.node_eccen > eccThreshold)];

% Which electrodes to plot? (Each electrode gets a subplot)
selectedElectrodes = setdiff(ElecsToInclude,ElecsToExclude);

% Sort electrodes by visual region (low to high)
selectedElectrodes = ecog_sortElectrodesonVisualArea(selectedElectrodes,trials.viselec.benson14_varea);

% Select only those channels with HRF broadband response higher than a certain threshold
whichTrials = {'HRF'};
[out] = ecog_plotTimecourses(trials, selectedElectrodes, whichTrials); % 'out' contains the plotted time courses and SEs

maxBb = [];
for ee = 1:length(selectedElectrodes)
    maxBb(ee) = max(out.broadband.(selectedElectrodes{ee}).mn);
end

elecNames = selectedElectrodes(maxBb > bbThreshold);

whichTrials = {'HRF'};
[out] = ecog_plotTimecourses(trials, elecNames, whichTrials);

%% Write out reduced version of trials for selected electrodes only

elecInx = ecog_matchChannels(elecNames, trials);

trials.channels  = trials.channels(elecInx,:);
trials.broadband = trials.broadband(elecInx,:,:);
trials.evoked    = trials.evoked(elecInx,:,:);

writeDir = '/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess';

saveName = fullfile(writeDir, sprintf('umcu%s_preproc_selectelecs',sub_label));
save(saveName, 'trials'); 

