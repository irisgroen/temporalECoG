
%% Get data
dataDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/preprocessed';

sub_label = {'som648'; 
             'som661'; 
             'som674'; 
             'som692'; 
             'som708'; 
             'som718';
             'som723';
             'umcuchaam';
             'umcubeilen'};

ses_label = {{'nyuecog01'}; 
             {'nyuecog01'};
             {'nyuecog01'}; 
             {'nyuecog01'}; 
             {'nyuecog01','nyuecog02'};
             {'nyuecog01'}; 
             {'nyuecog01', 'nyuecog02'}; 
             {'UMCUECOGday03'}; 
             {'UMCUECOGday03', 'UMCUECOGday06'}};
a = [];

for k = 1 : length(sub_label)
    for j = 1 : length(ses_label{k})
        D = bairanalysis_getepochs_visualelecs(dataDir, sub_label{k}, ses_label{k}{j});
        if j == 1
            a{k} = D;
        else
            % Check channels, concatenate sessions
            chan_names = D.trials.channels.name;
            assert(isequal(a{k}.trials.channels.name, chan_names))

            if isfield(summary(D.trials.channels), 'status')
                chan_status = D.trials.channels.status;
                assert(isequal(a{k}.trials.channels.status, chan_status))
                % todo: if this assertion fails, meaning that a visual
                % electrode is bad in one session but not another, include
                % only the good session, e.g. put one set of trials to nans
            end
                        
            a{k}.trials.broadband = cat(3,a{k}.trials.broadband,D.trials.broadband);
            a{k}.trials.evoked = cat(3,a{k}.trials.evoked,D.trials.evoked);
            a{k}.trials.events = [a{k}.trials.events; D.trials.events];
            
            a{k}.blank_trials.broadband = cat(3,a{k}.blank_trials.broadband,D.blank_trials.broadband);
            a{k}.blank_trials.evoked = cat(3,a{k}.blank_trials.evoked,D.blank_trials.evoked);
            a{k}.blank_trials.events = [a{k}.blank_trials.events; D.blank_trials.events];
            
            a{k}.ses = [a{k}.ses '+' D.ses];
        end
    end
end

% Exclude prf trials (not good for baseline)
for k = 1 : length(a)
    trialIdx = [];
    if max(contains(a{k}.trials.events.Properties.VariableNames, 'task_name'))
        trialIdx = ~contains(a{k}.trials.events.task_name, 'prf');
    else
        trialIdx = ~contains(a{k}.trials.events.trial_name, {'PRF', 'BLANK'});
    end
    
    a{k}.trials.events = a{k}.trials.events(trialIdx,:);
    a{k}.trials.broadband= a{k}.trials.broadband(:,:,trialIdx);
    a{k}.trials.evoked= a{k}.trials.evoked(:,:,trialIdx);
end

saveName = 'preproc_epoched_visualelecs.mat';
save(fullfile('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess', saveName), 'a', '-v7.3');

