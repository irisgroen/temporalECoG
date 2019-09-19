
function [out] = tde_getData(loadPrecomputed, inputDir, outputDir, subjectList, sessionList, epochTime)

% Description: 
%
% function [out] = tde_getData(loadPrecomputed, [inputDir], [outputDir], [subjectList], [sessionList], [epochTime])
%
% Input (use varargin?)
% - inputDir
% - outputDir
% - subjectslist
% - sessionList
% - epochTime
%
% Output
% a cell array with for each cell a struct with the following fields:
% - name
% - epochs
% - channels
% - events
% - t
%

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
    
    % Determine if we're loading or computing the data
    
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
        
        %% STEP 4: Save out a single preproc file for each subject 
        
        % Remove irrelevant columns from channels and events tables 
        events = removevars(events,{'onset','event_sample'});
        channels = removevars(channels,'notch');
        if isfield(summary(events),'stim_file'), removevars(events,'stim_file');end
        if isfield(summary(channels),'description'), removevars(channels,'description');end
        
        % Add a subject index column to both:
        events.subject_name = repmat({subject}, [height(events),1]);
        channels.subject_name = repmat({subject}, [height(channels),1]);
        
        % Save out the data
        fprintf('[%s] Saving data for subject %s to %s \n',mfilename, subject, outputDir);
        save(fullfile(outputDir, sprintf('%s_data_visualelecs.mat', subject)), 'epochs', 'channels', 'events', 't')
        
        % Collect into an output struct
        out{ii}.epochs   = epochs;
        out{ii}.channels = channels;
        out{ii}.events   = events;
        out{ii}.t        = t;
        
    end   
end

end





       
