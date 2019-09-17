
function [epochs, channels, events, t] = tde_getData(loadPrecomputed, inputDir, outputDir, subjectList, sessionList, epochLength)

% Description: 
%
% function [epochs, channels, events, t] = tde_getData(loadPrecomputed, [inputDir], [outputDir], [subjectList], [sessionList], [epochLength])
%
% Input
% - inputDir
% - outputDir
% - subjectslist? 
% - sessionList?
% - epochLength
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
    subjectList = {'som648'; 
                   'som661'; 
                   'som674'; 
                   'som692'; 
                   'som708'; 
                   'som718';
                   'som723';
                   'umcuchaam';
                   'umcubeilen'};
end

% <sessionList>
if ~exist('sessionList', 'var') || isempty(sessionList)
    sessionList = []; % default: all sessions
end

% <epochLength>
if ~exist('epochLength', 'var') || isempty(epochLength)
    epochLength = [-0.5 2];
end

% <tasks> (fixed for TDE project)
tasks = {'spatialpattern', 'temporalpattern', 'prf', 'soc'};

% <preprocessing data type> (fixed for TDE project)
description = 'broadband';

%% Loop across subjects

for ii = 1 : length(subjectList)
    
    subject = subjectList{ii};

    if loadPrecomputed
        
        %load from outputDir
        [epochs, channels, events, t] = load(fullfile(outputDir, sprintf('%s_data.mat', subject)));
        fprintf('[%s] Loading data from subject \n',mfilename, subject);
    
    else
        
        fprintf('[%s] Computing data from subject %s \n',mfilename, subject);

        % STEP 0: GET visual area matches for this subject
        fprintf('[%s] Computing matches with visual atlases...\n',mfilename);
        specs = [];
        specs.pID           = subject; 
        specs.plotmesh      = 'none';
        visualelectrodes    = electrode_to_nearest_node(specs);
 
        % STEP 1: GET THE DATA

        [data, channels, events, json] = bidsEcogGetPreprocData(inputDir, subjectList{ii}, sessionList, tasks, [], description);

        % STEP 2: SELECT A SUBSET OF CHANNELS 

        % Add visual area names (W and B) ecc, angle, sigma to channels table
        [channels] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes);

        % Make selection on visual only, index into data + channels
        chan_idx1 = find(~contains(channels.wangarea, 'none'));
        chan_idx2 = find(~contains(channels.bensonarea, 'none'));
        chan_idx = unique([chan_idx1; chan_idx2]); 
        
        data = data(chan_idx,:);
        channels = channels(chan_idx,:);
        
        % STEP 3: EPOCH THE DATA      
        [epochs, t] = ecog_makeEpochs(data, events.event_sample, epochLength, json.SamplingFrequency);  
        
        % Checks (Temporary)
%         trials = [];
%         trials.broadband = epochs;
%         trials.events = events;
%         trials.bb_bands = [50 200];
%         trials.time = t;
%         trials.channels = channels;
%         trials.viselec = visualelectrodes;
%         whichElectrodes = {'MO01', 'MO02', 'MO03', 'MO04'};
%         trialType = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%         %trialType = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%         specs.baselineType = 'selectedtrials';
%         specs.dataTypes = {'broadband'};
%         ecog_plotTimecourses(trials, whichElectrodes, trialType, specs)
%         
        % STEP 4 COMPUTE A BASELINE CORRECTION per trial per run
        
        % STEP 5 CHECK FOR OUTLIERS EPOCHS/CHANNELS
        % compute sum over time for each epoch?
        
    end

    % Save out a single preproc file for each subject that has
        % - epoched data (multiple sessions)
        % - reduced channel table, but with visual area matches
        % - events table (reduced?)
      
    % save(fullfile(outputDir, sprintf('%s_data_visual.mat', subject)), 'epochs', 'channels', 'events')
    
    % Concatenate epochs, channels, events across subjects?
    
    % Add subject index to event table
    
end





       
