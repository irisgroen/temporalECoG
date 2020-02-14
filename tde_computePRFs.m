function tde_computePRFs(recomputeData, doPlots)

% [results] = tde_computePRFs(recomputeData, recomputeFits, doPlots) 

% <saveDir> path to save parameters and fits; if empty, results are not
%   saved (default)
% <saveName> string to add to the save filename, if results are saved
%   (default empty)

% <recomputeData>
if ~exist('recomputeData','var') || isempty(recomputeData)
    recomputeData = false; % boolean
end

% <doPlots>
if ~exist('doPlots','var') || isempty(doPlots)
    doPlots = false; % boolean
end

if ~exist('saveDir', 'var'), saveDir = fullfile(analysisRootPath, 'results'); end
if ~exist('resultsStr', 'var'), resultsStr = 'prfs'; end

%% Get the PRF data
recomputeFlag = recomputeData;
subjects      = []; % will default to all subjects in subjectList.tsv
sessions      = []; % will default to all sessions per subject
tasks         = {'prf'};
description   = []; % will default to broadband
epochTime     = [-0.2 0.6];
sampleRate    = []; % will default to 512
saveStr       = 'prfdata';

[fulldata] = tde_getData(recomputeFlag, subjects, sessions, tasks, description, epochTime, sampleRate, saveStr);

% Load the stimulus apertures
stimName = fullfile(tdeRootPath, 'prf_apertures', 'bar_apertures.mat');
load(stimName, 'bar_apertures');
bar_apertures = imresize(bar_apertures, [100 100], 'nearest');

% Define a time window over which to average the broadband timecourse
time_win  = [0.05 0.55];

% Compute PRF timecourses for each subject
nSubjects = length(fulldata);

data2fit = cell(nSubjects,1);
 
for ii = 1:nSubjects
    
    subject = fulldata{ii}.subject;
    
    % Fix issues related to stimulus coding in some of the older datasets
    switch subject
        case 'beilen'
            % in this subject, trigger 41 was not sent because it was the
            % same as used for the blank (stimulus coding error).
            stimInx = setdiff(1:224, 41);
        case 'som661'
            % in this subject, more triggers were sent than the actual prf
            % bar positions --> need to interpolate? skip for now
            continue
        otherwise
            stimInx = 1:224;
    end
    
    % Compute average broadband response in time window
    t      = fulldata{ii}.t;
    epochs = fulldata{ii}.epochs;
    
    trials = squeeze(mean(epochs(t>time_win(1) & t<time_win(2),:,:),1)); 

    % Transpose to have channels in first dimension
    trials = trials';
    
    % Reshape to separate individual runs
    
    [nChans, nTrials] = size(trials);
    nStim = length(stimInx);
    
    nRuns = nTrials/nStim;
    
    ts = reshape(trials,[nChans nStim nRuns]);
    data2fit{ii}.ts = ts;
    data2fit{ii}.stim_inx = stimInx;
    
    % Make plots of the trials and of the PRF timecourses
    if doPlots
        % TO DO
    end
    
end
    
% Fit the PRF time courses with analyzePRF
tr             = 1;
opt.hrf        = 1;
opt.maxpolydeg = 0;
opt.xvalmode   = 0; 
opt.display    = 'off';

% Loop over subjects, fit data
for ii = 1:nSubjects

    subject = fulldata{ii}.subject;
    data = data2fit{ii};
    
    if ~isempty(data)
        
        % Average runs
        data = mean(data.ts,3);
    
        % Define stimulus
        stim_inx = data2fit{ii}.stim_inx;
        stimulus = bar_apertures(:,:,stim_inx);
        %stimulus = {bar_apertures,bar_apertures,bar_apertures,bar_apertures};

        results = analyzePRF(stimulus,data,tr,opt);

        % Save fits to results directory
        if ~isempty(saveDir)

            if ~exist(saveDir, 'dir'); mkdir(saveDir); end

            saveName = sprintf('%s_%s', subject, resultsStr);
            saveName = fullfile(saveDir, saveName);
            fprintf('[%s] Saving results to %s \n', mfilename, saveName);

            if exist(sprintf('%s.mat',saveName),'file')
                warning('[%s] Results file already exists! Writing new file with date-time stamp.',mfilename);
                saveName = sprintf('%s_%s', saveName, datestr(now,30));
                fprintf('[%s] Saving results to %s \n', mfilename, saveName);
            end
            save(saveName, 'data','results','tr','opt','subject');  
        end
    
        % Make plots of the estimated PRFs and PRF fits

        if doPlots
            % TO DO
        end
    end
end