function [data] = tde_getData(compute, inputDir, outputDir, subjectList, sessionList, epochTime, sampleRate)

% Description: 
%
% function [data] = tde_getData(compute, [inputDir], [outputDir], [subjectList], [sessionList], [epochTime], [sampleRate])
%
% Input (use varargin?)
% - inputDir
% - outputDir
% - subjectslist
% - sessionList
% - epochTime
% - sampleRate
%
% Output
% a cell array with for each cell a struct with the following fields:
% - name
% - epochs
% - channels
% - events
% - t

%% Define inputs 

% <compute>
if ~exist('compute', 'var') || isempty(compute)
	error('Please specify whether to compute the data (1) or to load from disk (0)');
end  

% <inputDir>
if ~exist('inputDir', 'var') || isempty(inputDir)
    %inputDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGBroadband';
    inputDir = fullfile(bidsRootPath, 'derivatives', 'ECoGBroadband');
end  

% <outputDir>
if ~exist('outputDir', 'var') || isempty(outputDir)
    %outputDir = '/Volumes/server/Projects/BAIR/Papers/TemporalDynamicsECoG/data';
	outputDir = fullfile(analysisRootPath, 'data');
end 

% <subjectList>
if ~exist('subjectList', 'var') || isempty(subjectList)
    subjectList_fname = fullfile(tdeRootPath, 'subjectlist.tsv');
    T = readtable(subjectList_fname, 'FileType', 'text');
    subjectList = T.participant_id;
end

% <sessionList>
if ~exist('sessionList', 'var') || isempty(sessionList)
    sessionList = []; % default: all sessions
end

% <epochTime>
if ~exist('epochTime', 'var') || isempty(epochTime)
    epochTime = [-0.2 1];
end

% <sampleRate>
if ~exist('sampleRate', 'var') || isempty(sampleRate)
    sampleRate = 512;
end

% <tasks> (fixed for TDE project)
tasks = {'spatialpattern', 'temporalpattern', 'soc'};

% <preprocessing data type> (fixed for TDE project)
description = 'broadband';

%% Loop across subjects
data = cell(length(subjectList),1);

for ii = 1 : length(subjectList)
    
    subject = subjectList{ii};
    
    % Determine if we're loading or computing the data
    
    if ~compute
        
        % load from outputDir    
        data{ii} = load(fullfile(outputDir, sprintf('%s_data_visualelecs.mat', subject)));
        fprintf('[%s] Loading data for subject %s \n',mfilename, subject);
    
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

        [subdata, channels, events] = bidsEcogGetPreprocData(inputDir, subject, sessionList, tasks, [], description);

        %% STEP 2: SELECT A SUBSET OF CHANNELS 

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
        fprintf('[%s] Step 5: Saving data for subject %s to %s \n',mfilename, subject, outputDir);
        save(fullfile(outputDir, sprintf('%s_data_visualelecs.mat', subject)),'subject', 'epochs', 't', 'events', 'channels')
        
        % Collect into an output struct
        data{ii}.subject  = subject;
        data{ii}.epochs   = epochs;
        data{ii}.t        = t;
        data{ii}.events   = events;
        data{ii}.channels = channels;

    end   
end
fprintf('[%s] Done! \n',mfilename);
end





       
