
%% function [epochs, channels, events] = tde_getData(recomputeflag 0/1)

% function [epochs, channels, events] = tde_getData(recomputeflag 0/1, [dataDir], [subjectList], [sessionList], [epochLength], [inputDataType])

% make inputs as well??
% - subjectslist? sessionList?
% - epochLength
% - inputdataType (broadband)

dataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoGBroadband';

sub_label = {'som648'; 
             'som661'; 
             'som674'; 
             'som692'; 
             'som708'; 
             'som718';
             'som723';
             'umcuchaam';
             'umcubeilen'};

% ses_label = {{'nyuecog01'}; 
%              {'nyuecog01'};
%              {'nyuecog01'}; 
%              {'nyuecog01'}; 
%              {'nyuecog01','nyuecog02'};
%              {'nyuecog01'}; 
%              {'nyuecog01','nyuecog02'}; 
%              {'UMCUECOGday03'}; 
%              {'UMCUECOGday03','UMCUECOGday06'}};

tasks = {'spatialpattern', 'temporalpattern', 'prf'};
epochlength = [-0.5 2];


for k = 1 : length(sub_label)
    
     
    % GET visual matches for this subject
        
    % viselec = electrode_to_nearest_node
    
    % bidsSpecifySessions
    
    % STEP 1: GET THE DATA

    % [data, events, chans] = bidsEcogGetPreprocData(dataPath, sub_label{k}, [], tasks)
   
    % STEP 2: SELECT A SUBSET OF CHANNELS 

    % only those with matches to visual regions

    % [channels] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes)
    % make selection, index into data + channels

    % STEP 3: EPOCH THE DATA

    % should segment the trials according to specific epoch length.
    % input: event onsets (samples?) 

    % [epocheddata, events, chans] = ecog_makeEpochs()

    % Save out a single preproc file for each subject (containing both
    % sessions) that has
    % - epoched data
    % - reduced channel table, but with visual area matches
    % - events table or just list of trial numbers/names?
    
end





       
