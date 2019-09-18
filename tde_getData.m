
function [out] = tde_getData(loadPrecomputed, inputDir, outputDir, subjectList, sessionList, epochTime, baselineTime)

% Description: 
%
% function [epochs, channels, events, t] = tde_getData(loadPrecomputed, [inputDir], [outputDir], [subjectList], [sessionList], [epochTime], [baselineTime])
%
% Input (use varargin?)
% - inputDir
% - outputDir
% - subjectslist? 
% - sessionList?
% - epochTime
% - baselineTime
%
% Output
%
% Example

%% Define inputs 

% <compute>
if ~exist('loadPrecomputed', 'var') || isempty(loadPrecomputed)
	error('Please specify whether to load the precomputed data (1) or to compute the data anew from the BIDS directory (0)');
end  

% <inputDir>
if ~exist('inputDir', 'var') || isempty(inputDir)
    inputDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGBroadband';
end  

% <outputDir>
if ~exist('ouputDir', 'var') || isempty(outputDir)
    outputDir = '/Volumes/server/Projects/BAIR/Papers/TemporalDynamicsECoG/data';
end 

% <subjectList>
if ~exist('subjectList', 'var') || isempty(subjectList)
    subjectList = {'chaam';
                   'beilen';
                   'som648'; 
                   'som661'; 
                   'som674'; 
                   'som692'; 
                   'som708'; 
                   'som718';
                   'som723';
                   };
end

% <sessionList>
if ~exist('sessionList', 'var') || isempty(sessionList)
    sessionList = []; % default: all sessions
end

% <epochTime>
if ~exist('epochTime', 'var') || isempty(epochTime)
    epochTime = [-0.5 2];
end

% <baselineTime>
if ~exist('baselineTime', 'var') || isempty(baselineTime)
    baselineTime = [-0.2 0];
end

% <tasks> (fixed for TDE project)
tasks = {'spatialpattern', 'temporalpattern', 'bairspatialpattern', 'bairtemporalpattern','soc'};

% <preprocessing data type> (fixed for TDE project)
description = 'broadband';

%% Loop across subjects

for ii = 1 : length(subjectList)
    
    subject = subjectList{ii};
    out{ii}.name = subject;
    
    if loadPrecomputed
        
        % load from outputDir    
        out{ii} = load(fullfile(outputDir, sprintf('%s_data_visualelecs.mat', subject)));
        fprintf('[%s] Loading data for subject %s \n',mfilename, subject);
    
    else
        
        fprintf('[%s] Computing data for subject %s \n',mfilename, subject);

        % STEP 0: GET visual area matches for this subject
        fprintf('[%s] Computing matches with visual atlases...\n',mfilename);
        specs = [];
        specs.pID           = subject; 
        specs.plotmesh      = 'none';
        BIDSformatted       = 1;
        visualelectrodes    = electrode_to_nearest_node(specs, BIDSformatted);
 
        %% STEP 1: GET THE DATA

        [data, channels, events] = bidsEcogGetPreprocData(inputDir, subjectList{ii}, sessionList, tasks, [], description);

        %% STEP 2: SELECT A SUBSET OF CHANNELS 

        % Add visual area names (W and B) ecc, angle, sigma to channels table
        [channels] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes);

        % Make selection on visual only, index into data + channels
        fprintf('[%s] Selecting channels with visual matches...\n',mfilename);
        
        %chan_idx1 = find(contains(channels.type,'ecog') & ~contains(channels.wangarea, 'none'));
        %chan_idx2 = find(contains(channels.type,'ecog') & ~contains(channels.bensonarea, 'none'));
        chan_idx1 = find(~contains(channels.wangarea, 'none'));
        chan_idx2 = find(~contains(channels.bensonarea, 'none'));
        
        chan_idx = unique([chan_idx1; chan_idx2]);
        data = data(chan_idx,:);
        channels = channels(chan_idx,:);
        
        fprintf('[%s] Found %d channels with visual matches out of %d ecog channels \n', ...
            mfilename, length(chan_idx), length(find(contains(lower(channels.type), {'ecog', 'seeg'}))));
                
        %% STEP 3: EPOCH THE DATA 
        fprintf('[%s] Epoching data...\n',mfilename);
        
        [epochs, t] = ecog_makeEpochs(data, events.event_sample, epochTime, channels.sampling_frequency(1));  
        
        fprintf('[%s] Found %d epochs across %d runs and %d sessions \n', ...
            mfilename, size(epochs,3), length(unique(events.run_name)), length(unique(events.session_name)));
       
        %% STEP 4 CHECK FOR OUTLIERS EPOCHS/CHANNELS
        % fprintf('[%s] Checking for bad epochs ...\n',mfilename);

        % compute sum over time for each epoch?
        % remove epochs that have 20x the average
        % if more than X epochs for a channel, remove entire channel
        % output a description of how many trials were removed (write to
        % file?)
        
        %% STEP 5 convert to percent signal change
        fprintf('[%s] Converting epochs to percent signal change...\n',mfilename);

        % Provide run index to perform separately for each run and session
        [~,~,ses_idx]= unique(events.session_name);
        [~,~,run_idx] = unique(events.run_name);
        idx = (ses_idx*100)+run_idx;
        
        [epochs] = ecog_normalizeEpochs(epochs, t, baselineTime, 'percentsignalchange', idx);
        
        % Remove onsets from events (does not apply to epochs)
        events = removevars(events,{'onset','event_sample'});
    
        % Checks (Temporary; TO write simpler function e.g. bidsEcogPlotEpochs)
%         trials = [];
%         trials.broadband = epochs;
%         trials.events = events;
%         trials.bb_bands = [50 200];
%         trials.time = t;
%         trials.channels = channels;
%         trials.viselec = visualelectrodes;
%         whichElectrodes = channels.name;
%         specs.baselineType = 'selectedtrials';
%         specs.dataTypes = {'broadband'};
%         trialType = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%         ecog_plotTimecourses(trials, whichElectrodes, trialType, specs);
%         trialType = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%         ecog_plotTimecourses(trials, whichElectrodes, trialType, specs);
%         trialType = {'CRF-1', 'CRF-2', 'CRF-3', 'CRF-4', 'CRF-5'};
%         ecog_plotTimecourses(trials, whichElectrodes, trialType, specs);
        
        %% STEP 6 Save out a single preproc file for each subject that has
        fprintf('[%s] Saving data for subject %s to %s \n',mfilename, subject, outputDir);
        save(fullfile(outputDir, sprintf('%s_data_visualelecs.mat', subject)), 'epochs', 'channels', 'events', 't')
            
        out{ii}.epochs   = epochs;
        out{ii}.channels = channels;
        out{ii}.events   = events;
        out{ii}.t        = t;
        
    end   
end

end





       
