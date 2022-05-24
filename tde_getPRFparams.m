function channels = tde_getPRFparams(channels, loadDir, resultsStr)
% Loads prf parameters for each electrode in channels table as saved out 
% by tde_fitPRFs.
%
% channels = tde_getPRFparams(channels, [loadDir], [resultsStr])
%
% <loadDir> path to saved parameters and fits
%   default: fullfile(analysisRootPath, 'prfs')
% <resultsStr> string added to the prf filename 
%   default: 'prfs'
%
% see tde_fitPRFs.m
%
% 2020 Iris Groen

% <loadDir>
if ~exist('loadDir', 'var') || isempty(loadDir), loadDir = fullfile(analysisRootPath, 'prfs'); end

% <resultsStr>
if ~exist('resultsStr', 'var') || isempty(resultsStr), resultsStr = 'prffits'; end

if ~isfield(summary(channels),'subject_name')
    error('No subject info in channels table, cannot look up corresponding PRF data')
end

subjectList = unique(channels.subject_name);
nSubjects = length(subjectList);

channels_with_prf = [];
sortorder = nan(height(channels),1);

for ii = 1:nSubjects
    
    subject = subjectList{ii};
    chan_idx = contains(channels.subject_name, subject);  
    sortorder(chan_idx) = ii;
    
    loadName = sprintf('sub-%s_%s.mat', subject, resultsStr);
	loadName = fullfile(loadDir, loadName);
    if ~isfile(loadName)
        
        fprintf('[%s] Could not locate prf results for subject %s using %s \n', mfilename, subject, loadName);
        [channels_new] = bair_addAnalyzePrfResultsToChannelTable(channels(chan_idx,:),[]);

    else
    
        fprintf('[%s] Loading prf results from %s \n', mfilename, loadName);
        prf = load(loadName);        
        [channels_new] = bair_addAnalyzePrfResultsToChannelTable(channels(chan_idx,:),prf.results);
        
    end
    
    % concatenate
    channels_with_prf = vertcat(channels_with_prf, channels_new);
   
end

channels = sortrows(channels_with_prf, sortorder);

end