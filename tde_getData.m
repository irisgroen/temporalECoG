function [data] = tde_getData(compute, subjects, sessions, tasks, description, epochTime, sampleRate, saveStr, saveDir, dataDir )

% Read in ECoG data from visual template matching channels for multiple
% subjects for a given set of tasks.
%
% [data] = tde_getData(compute, [subjects], [sessions], [tasks], ...
%                      [description], [epochTime], [sampleRate], ...
%                      [saveStr], [saveDir], [dataDir]
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
% - description : string to match to the bids field desc-<desc>
%          in the input filenames. default: 'broadband';
% - epochTime : [t_start t_stop] array defining the epoch window
%          default: [-0.2 1];
% - sampleRate : desired sample rate in Hz for all datasets.
%          Datasets with rates will be downsampled. default: 512            
% - saveStr : string to be added to filename for saved out data
%          default: 'tdedata'
% - saveDir : directory to write data to 
%          default: fullfile(analysisRootPath, 'data');
% - dataDir : directory to get data from 
%          default: fullfile(bidsRootPath, 'derivatives', 'ECoGBroadband')
%
% OUTPUT
% A cell array with for each cell a struct with the following fields:
% - name
% - epochs
% - channels
% - events
% - t
%
% NOTES
% - Data should be bids-formatted.
% - Function will perform the following steps:
%   STEP 0: Match electrode positions to wang and benson atlases
%   STEP 1: Read in the datafiles
%   STEP 2: Select channels with a visual match to either of the atlases.
%   STEP 3: Downsample data (if higher sample rate than SampleRate) and
%           shift onsets for UMCU patients.
%   STEP 4: Epoch the data according to the onsets in the events.tsv files
%           found in the dataDir according to epochTime
%   STEP 5: Save out data for each subject in the saveDir.
% 
% Uses electrode_to_nearest_node.m bidsEcogGetPreprocData.m
%       ecog_makeEpochs.m bair_addVisualAtlasNamesToChannelTable
% 
% IG 2020

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

% <description> 
if ~exist('description', 'var') || isempty(description)
    description = 'broadband'; % TDE
end

% <epochTime>
if ~exist('epochTime', 'var') || isempty(epochTime)
    epochTime = [-0.2 1];
end

% <sampleRate>
if ~exist('sampleRate', 'var') || isempty(sampleRate)
    sampleRate = 512;
end

% <saveStr>
if ~exist('saveStr', 'var') || isempty(saveStr)
	saveStr = 'tdedata';
end 

% <saveDir>
if ~exist('saveDir', 'var') || isempty(saveDir)
	saveDir = fullfile(analysisRootPath, 'data');
end 

% <dataDir>
if ~exist('dataDir', 'var') || isempty(dataDir)
    dataDir = fullfile(bidsRootPath, 'derivatives', 'ECoGBroadband');
end 

%% Loop across subjects
if ~iscell(subjects), subjects = {subjects}; end

data = cell(length(subjects),1);

for ii = 1 : length(subjects)
    
    subject = subjects{ii};
    
    % Determine if we're loading or computing the data
    
    if ~compute
        
        % load from outputDir
        fileName = fullfile(saveDir, sprintf('%s_%s_visualelecs.mat', subject, saveStr));
        if exist(fileName, 'file')
            data{ii} = load(fileName);
            fprintf('[%s] Loading data for subject %s \n',mfilename, subject);
        end
    
    else
        
        fprintf('[%s] Computing data for subject %s \n',mfilename, subject);

        %% STEP 0: GET visual area matches for this subject
        fprintf('[%s] Step 0: Computing matches with visual atlases...\n',mfilename);
        
        specs = [];
        specs.pID           = subject; 
        specs.plotmesh      = 'none';
        BIDSformatted       = 1;
        visualelectrodes    = electrode_to_nearest_node(specs, BIDSformatted);
 
        %% STEP 1: GET THE DATA
        fprintf('[%s] Step 1: Loading data...\n',mfilename);
        
        if ~isempty(sessions), session = sessions{ii}; else, session = []; end

        [subdata, channels, events] = bidsEcogGetPreprocData(dataDir, subject, session, tasks, [], description);

        if isempty(subdata)
            
            % We did not find any data for this subject; move to next one.
            warning('No data found for subject %s!', subject);
            continue
            
        %% STEP 2: SELECT A SUBSET OF CHANNELS 
        
        
        else
            % Add visual area names (W and B) ecc, angle, sigma to channels table
            [channels] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes);

            % Make selection on visual only, index into data + channels
            fprintf('[%s] Step 2: Selecting channels with visual matches \n',mfilename);

            chan_idx1 = find(~contains(channels.wangarea, 'none') & contains(channels.status, 'good'));
            chan_idx2 = find(~contains(channels.bensonarea, 'none') & contains(channels.status, 'good'));        
            chan_idx = unique([chan_idx1; chan_idx2]);

            fprintf('[%s] Step 2: Found %d channels with visual matches out of %d ecog channels \n', ...
                mfilename, length(chan_idx), length(find(contains(lower(channels.type), {'ecog', 'seeg'}))));

            subdata = subdata(chan_idx,:);
            channels = channels(chan_idx,:);

            %% STEP 3: CHECK SAMPLE RATES ANDS SHIFT DATA

            fprintf('[%s] Step 3: Checking sample rates...\n',mfilename);

            if channels.sampling_frequency(1) ~= sampleRate
                fprintf('[%s] Step 3: Sample rate does not match requested sample rate. Resampling \n',mfilename);
                subdata = downsample(subdata', channels.sampling_frequency(1)/sampleRate)';
                events.event_sample = round(events.event_sample/(channels.sampling_frequency(1)/sampleRate));
                channels.sampling_frequency(:) = sampleRate;
            end

            % SHIFT the UMCU data 
            if contains(subject, {'chaam', 'beilen'}) % SHOULD BE READ FROM participants.tsv, if site column = umcu
                fprintf('[%s] Step 3: This is a umcu patient. Shifting onsets \n',mfilename);

                % Shift onsets
                shiftInSeconds = 0.072; % 72 ms; determined through cross correlation, see s_determineOnsetShiftUMCUvsNYU.m
                shiftInSamples = round(shiftInSeconds*channels.sampling_frequency(1)); 
                events.onset = events.onset + shiftInSeconds;
                events.event_sample = events.event_sample + shiftInSamples; 

                % Remove electrodes identified as epileptic % TEMPORARY UNTIL GIO
                % FIXES BIDS FORMATTING FOR THIS PATIENT AFTER WHICH THESE
                % ELECTRODES SHOULD BE INDICATED AS "BAD" IN THE ELECTRODES FILES
                if contains(subject,'chaam')
                    chan_idx = find(~contains(channels.name, {'Oc12', 'Oc13', 'Oc21', 'Oc22'}));
                    subdata = subdata(chan_idx,:);
                    channels = channels(chan_idx,:);
                end
            end

            %% STEP 4: EPOCH THE DATA 

            fprintf('[%s] Step 4: Epoching data \n',mfilename);

            [epochs, t] = ecog_makeEpochs(subdata, events.event_sample, epochTime, channels.sampling_frequency(1));  

            fprintf('[%s] Step 4: Found %d epochs across %d runs and %d sessions \n', ...
                mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));

            %% STEP 5: Save out a single preproc file for each subject 

            % Remove irrelevant/redundant columns from channels and events tables 
            if isfield(summary(events),'onset'), events = removevars(events,'onset');end
            if isfield(summary(events),'event_sample'), events = removevars(events,'event_sample');end
            if isfield(summary(events),'stim_file'), events = removevars(events,'stim_file');end
            if isfield(summary(channels),'notch'), channels = removevars(channels,'notch');end
            if isfield(summary(channels),'status'), channels = removevars(channels,'status');end
            if isfield(summary(channels),'description'), channels = removevars(channels,'description');end
            if isfield(summary(channels),'status_description'), channels = removevars(channels,'status_description');end

            % Add a subject index column to channels and events tables:
            events.subject_name = repmat({subject}, [height(events),1]);
            channels.subject_name = repmat({subject}, [height(channels),1]);

            % Save out the data
            fprintf('[%s] Step 5: Saving data for subject %s to %s \n',mfilename, subject, saveDir);
            saveName = sprintf('%s_%s_visualelecs.mat', subject, saveStr);
            saveName = fullfile(saveDir, saveName);
            save(saveName,'subject', 'epochs', 't', 'events', 'channels')

            % Collect into an output struct
            data{ii}.subject  = subject;
            data{ii}.epochs   = epochs;
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

fprintf('[%s] Done! \n',mfilename);
end





       
