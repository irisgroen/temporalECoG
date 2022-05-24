function [data] = tde_getData(compute, subjects, sessions, tasks, epochTime, sampleRate, bidsDir, saveStr, saveDir)

% Read in ECoG voltage and broadband data from corresponding BIDS
% derivatives directories for specified subjects, sessions and tasks, and
% saves these out as a separate <sub>_data .mat files for each subject in
% <tdeRootPath>/analysis/data. If 'compute' is set to 0, data will not be
% re-extracted but instead loaded from this directory.
%
% When extracting data, the following operations are performed:
% - Select channels with matches to prespecified retinotopic atlases
% (Benson14 and Wang maximum probability atlas) will be retained.
% - Shift UMCU data by 72 ms
% - Epoch the data according to the onsets in events.tsv
% - Combine all sessions and tasks into a single data file per subject.
%
% [data] = tde_getData(compute, [subjects], [sessions], [tasks], [epochTime], ...
%                               [sampleRate], [bidsDir], [saveStr], [saveDir])
%
% INPUT (required):
% - compute : boolean indicating whether to compute or read from disk.
%
% INPUT (optional):
% - subjects : cell array of subject names. If empty, script
%           will read subject names from subject_id column in
%           subjectlist.tsv file in tdeRootPath.m
% - sessions : cell array of session names of the same 
%           dimensions as subjectList. default: all sessions 
% - tasks : cell array of task names to match to the bids field 
%          task-<taskname> in the input filenames.
%          default: {'spatialpattern', 'temporalpattern', 'soc'};
% - epochTime : [t_start t_stop] array defining the epoch window
%          default: [-0.2 1.2];
% - sampleRate : desired sample rate in Hz for all datasets.
%          Datasets with rates will be downsampled. default: 512            
% - bidsDir : directory to read data from
%          default: fullfile(bidsRootPath);
% - saveStr : string to be added to filename for saved out data
%          default: 'tdedata'
% - saveDir : directory to write data to 
%          default: fullfile(analysisRootPath, 'data');
%
% OUTPUT
% A cell array with for each cell a struct with the following fields:
% - subject  (string with subjectname)
% - epochs_v (time x events x channels matrix with voltage data)
% - epochs_b (time x events x channels matrix with broadband data)
% - channels (bids-formatted channel table)
% - events   (bids-formatted events table)
% - t        (vector with time points)
%
% NOTES
% - Data should be bids-formatted.
% - Function will perform the following steps:
%   STEP 0: Match electrode positions to wang and benson atlases
%   STEP 1: Read in the time series data: both broadband and voltage; also
%           resample if sample rate does not match SampleRate argument
%   STEP 2: Select channels with a visual match to either of the atlases.
%   STEP 3: Deal with UMCU patients: shift onsets.
%   STEP 4: Epoch the data according to the onsets in the events.tsv files
%           found in the dataDir according to epochTime
%   STEP 5: Save out data for each subject in the saveDir.
% 
% Uses electrode_to_nearest_node.m bidsEcogGetPreprocData.m
%       ecog_makeEpochs.m bair_addVisualAtlasNamesToChannelTable
% 
% 2020 Iris Groen

%% Define inputs 

% <compute>
if ~exist('compute', 'var') || isempty(compute)
	error('Please specify whether to compute the data (1) or to load from disk (0)');
end  

% <subjects>
if ~exist('subjects', 'var') || isempty(subjects)
    subjectList_fname = fullfile(tdeRootPath, 'subjectlist.tsv');
    T = readtable(subjectList_fname, 'FileType', 'text');
    subjects = T.participant_id;
end

% <sessions>
if ~exist('sessions', 'var') || isempty(sessions)
    sessions = []; % default: all sessions
end

% <tasks> 
if ~exist('tasks', 'var') || isempty(tasks)
    tasks = {'spatialpattern', 'temporalpattern', 'soc'}; % TDE 
end

% <epochTime>
if ~exist('epochTime', 'var') || isempty(epochTime)
    epochTime = [-0.1 1.2];
end

% <sampleRate>
if ~exist('sampleRate', 'var') || isempty(sampleRate)
    sampleRate = 512;
end

% <readDir>
if ~exist('bidsDir', 'var') || isempty(bidsDir)
	bidsDir = bidsRootPath;
end 

% <saveStr>
if ~exist('saveStr', 'var') || isempty(saveStr)
	saveStr = 'tdedata';
end 

% <saveDir>
if ~exist('saveDir', 'var') || isempty(saveDir)
	saveDir = fullfile(analysisRootPath, 'data');
end 

if ~exist(saveDir, 'dir'), mkdir(saveDir); end

%% Loop across subjects
if ~iscell(subjects), subjects = {subjects}; end

data = cell(length(subjects),1);

for ii = 1 : length(subjects)
    
    subject = subjects{ii};
    
    % Determine if we're loading or computing the data
    
    if ~compute
        
        % load from outputDir
        fileName = fullfile(saveDir, sprintf('sub-%s_%s.mat', subject, saveStr));
        if exist(fileName, 'file')
            data{ii} = load(fileName);
            fprintf('[%s] Loading data for subject %s \n',mfilename, subject);
        else
            fprintf('[%s] Could not locate datafile to load from disk for subject %s \n', mfilename,subject);
        end
        
    else
        
        fprintf('[%s] Computing data for subject %s \n',mfilename, subject);
  
        %% STEP 1: GET THE DATA
        fprintf('[%s] Step 1: Loading data...\n',mfilename);
        
        if ~isempty(sessions), session = sessions{ii}; else, session = []; end

        % Read in voltage data
        dataDir = fullfile(bidsDir, 'derivatives', 'ECoGCAR');
        [data_v, ~, ~] = bidsEcogGetPreprocData(dataDir, subject, session, tasks, [], 'reref', sampleRate);
        if isempty(data_v), warning('[%s] No voltage data found for subject %s!', mfilename, subject); end

        % Read in broadband data
        dataDir = fullfile(bidsDir, 'derivatives', 'ECoGBroadband');
        [data_b, channels, events] = bidsEcogGetPreprocData(dataDir, subject, session, tasks, [], 'broadband', sampleRate);
        if isempty(data_b), warning('[%s] No broadband data found for subject %s!', mfilename, subject); end
        
        if isempty(data_v) && isempty(data_b), continue; end
            
        % Read in electrode data and match to atlas
        dataDir = fullfile(bidsDir);
        atlasName = {'benson14_varea',  'wang15_mplbl', 'wang15_fplbl', 'benson14_eccen', 'benson14_angle', 'benson14_sigma'};
        [electrodes] = bidsEcogMatchElectrodesToAtlas(dataDir, subject, session, atlasName, [], 0);
        
        % Reduce data to electrodes with coordinates only
        chan_idx = ecog_matchChannels(electrodes.name, channels.name);        
        data_v = data_v(chan_idx,:);
        data_b = data_b(chan_idx,:);
        channels = channels(chan_idx,:);
        
        % Add electrode atlas info to channels
        assert(height(electrodes) == height(channels));
        assert(isequal(electrodes.name, channels.name));
        col_idx = ~contains(electrodes.Properties.VariableNames, channels.Properties.VariableNames);
        channels = [channels electrodes(:,col_idx)];
        
        %% STEP 2: SELECT A SUBSET OF CHANNELS 
           
        if ~isempty(data_v) && ~isempty(data_b)
            
            % Make selection on visual only, index into data + channels
            fprintf('[%s] Step 2: Selecting channels with visual matches \n',mfilename);
            
            chan_idx1 = find(~contains(channels.benson14_varea, 'none') & contains(channels.status, 'good'));
            chan_idx2 = find(~contains(channels.wang15_mplbl, 'none') & contains(channels.status, 'good'));        
            chan_idx = unique([chan_idx1; chan_idx2]);

            if ~isempty(chan_idx)
                fprintf('[%s] Step 2: Found %d channels with visual matches out of %d ecog channels \n', ...
                mfilename, length(chan_idx), length(find(contains(lower(channels.type), {'ecog', 'seeg'}))));
            else
                warning('No visual matches found for subject %s!\n', subject);
                continue
            end
            
            % Reduce data to selected channels only.
            data_v = data_v(chan_idx,:);
            data_b = data_b(chan_idx,:);
            channels = channels(chan_idx,:);
            
            %% STEP 3: DEAL WITH UMCU DATA (shift onsets)

            % SHIFT the UMCU data 
            if contains(subject, {'p01', 'p02'}) 
                fprintf('[%s] Step 3: This is a umcu patient. Shifting onsets \n',mfilename);

                % Shift onsets
                shiftInSeconds = 0.072; % 72 ms; determined through cross correlation, see s_determineOnsetShiftUMCUvsNYU.m
                events.onset = events.onset + shiftInSeconds;

            end

            %% STEP 4: EPOCH THE DATA 

            fprintf('[%s] Step 4: Epoching data \n',mfilename);

            [epochs_v, ~] = ecog_makeEpochs(data_v, events.onset, epochTime, channels.sampling_frequency(1));  
            [epochs_b, t] = ecog_makeEpochs(data_b, events.onset, epochTime, channels.sampling_frequency(1));  
            
            fprintf('[%s] Step 4: Found %d epochs across %d runs and %d sessions \n', ...
                mfilename, size(epochs_b,2), length(unique(events.run_name)), length(unique(events.session_name)));

            %% STEP 5: Save out a single preproc file for each subject 

            % Remove irrelevant/redundant columns from events table
            if isfield(summary(events),'onset'), events = removevars(events,'onset');end
            if isfield(summary(events),'stim_file'), events = removevars(events,'stim_file');end
            
            % Remove irrelevant/redundant columns from channels table
            if isfield(summary(channels),'notch'), channels = removevars(channels,'notch');end
            if isfield(summary(channels),'status'), channels = removevars(channels,'status');end
            if isfield(summary(channels),'description'), channels = removevars(channels,'description');end
            if isfield(summary(channels),'status_description'), channels = removevars(channels,'status_description');end
            if isfield(summary(channels),'size'), channels = removevars(channels,'size');end
            if isfield(summary(channels),'material'), channels = removevars(channels,'material');end
            if isfield(summary(channels),'manufacturer'), channels = removevars(channels,'manufacturer');end

            % Add a subject index column to channels and events tables:
            events.subject_name = repmat({subject}, [height(events),1]);
            channels.subject_name = repmat({subject}, [height(channels),1]);

            % Save out the data
            fprintf('[%s] Step 5: Saving data for subject %s to %s \n',mfilename, subject, saveDir);
            saveName = sprintf('sub-%s_%s.mat', subject, saveStr);
            saveName = fullfile(saveDir, saveName);
            save(saveName,'subject', 'epochs_b', 'epochs_v', 't', 'events', 'channels')

            % Collect into an output struct
            data{ii}.subject  = subject;
            data{ii}.epochs_b = epochs_b;
            data{ii}.epochs_v = epochs_v;
            data{ii}.t        = t;
            data{ii}.events   = events;
            data{ii}.channels = channels;
        end
    end   
end

% Remove empty cells from the output
emptycells=[];
for ii = 1:length(data)
    if isempty(data{ii}), emptycells = [emptycells ii]; end
end
data(emptycells) = [];

if ~isempty(data)
    fprintf('[%s] Done! \n',mfilename);
else
    fprintf('[%s] No data found! \n',mfilename);
end
end





       
