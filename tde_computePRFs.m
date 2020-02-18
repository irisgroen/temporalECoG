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

%% Set paths
if ~exist('saveDir', 'var'), saveDir = fullfile(analysisRootPath, 'prfs'); end
if ~exist('resultsStr', 'var'), resultsStr = 'prfs'; end
if doPlots
    plotSaveDir = fullfile(analysisRootPath, 'figures', 'prfs');
    if ~exist(plotSaveDir, 'dir'); mkdir(fullfile(plotSaveDir, 'data'));  mkdir(fullfile(plotSaveDir, 'modelfits'));end
end

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
    
    subject  = fulldata{ii}.subject;
    t        = fulldata{ii}.t;
    epochs   = fulldata{ii}.epochs;
    
    % Determine which stimuli to select for the PRF timecourse
    switch subject
        case 'beilen'
            % in this subject, trigger 41 was not sent because it was the
            % same as used for the blank (stimulus coding error).
            stimInx = setdiff(1:224, 41);
        case 'som661'
            % in this subject, more triggers were sent than the actual prf
            % bar positions --> need to interpolate? skip for now
            continue
        case 'som674'
            % the first two prf runs in this subject are bad (broke
            % fixation), we will not analyze these
            start_ind = (size(epochs,2)/2)+1;
            epochs = epochs(:,start_ind:end,:); stimInx = 1:224;
        otherwise
            stimInx = 1:224;
    end
    
    % Normalize epochs, within run
    [~, nTrials, nChans] = size(epochs);
    nStim = length(stimInx);    
    nRuns = nTrials/nStim; 
    run_indices = []; 
    for jj = 1:nRuns
        run_indices = [run_indices; ones(nStim,1)*jj];
    end
    %[epochs] = ecog_normalizeEpochs(epochs, t, [], [], run_indices);
    
    % Compute average broadband response in time window
    trials = squeeze(mean(epochs(t>time_win(1) & t<time_win(2),:,:),1)); 

    % Transpose to have channels in first dimension
    trials = trials';
    
    % Reshape to separate individual runs    
    ts = reshape(trials,[nChans nStim nRuns]);
    
    data2fit{ii}.ts       = ts;
    data2fit{ii}.stim_inx = stimInx;

    % Make plots of the trials and of the PRF timecourses
    if doPlots
        channels = fulldata{ii}.channels;
        events   = fulldata{ii}.events;
        nEpochs  = size(epochs,2);

        % Plot individual trials
        fprintf('[%s] Plotting prf trials for subject %s \n',mfilename, subject);
        figureName = sprintf('%s_prftrials', subject);
        figure('Name', figureName); plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
        for el = 1:nChans
            subplot(plotDim1,plotDim2,el); hold on
            for kk = 1:nEpochs
                %ecog_plotSingleTimeCourse(t, epochs(:,kk,el));
                plot(t, epochs(:,kk,el)); 
            end
            plot(t, squeeze(mean(epochs(:,kk,el),2)), 'k', 'LineWidth', 3); axis tight;

            yLim = get(gca, 'YLim');
            line([0 0], yLim,'LineStyle', ':', 'Color', 'k');
            line([t(1) t(end)], [0 0],'LineStyle', ':', 'Color', 'k');
            plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        
            title(plotTitle);
            axis tight
            if el == 1; xlabel('Time (s)'); ylabel('Broadband signal change');end
            %set(gcf, 'Position', [150 100 1500 1250]);
            set(gcf, 'Position', get(0,'screensize'));
            set(gca, 'FontSize', 14);
        end
        saveas(gcf, fullfile(plotSaveDir,'data', figureName), 'png'); close;
        
        % Plot PRF timecourses for each run + average
        fprintf('[%s] Plotting prf timecourses data for subject %s \n',mfilename, subject);
        figureName = sprintf('%s_prftimecourses', subject);
        figure('Name', figureName); plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
        for el = 1:nChans
            subplot(plotDim1,plotDim2,el); hold on
            plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        
            runNames = [];
            for kk = 1:nRuns
                plot(1:nStim, squeeze(ts(el,:,kk)), 'LineWidth', 1);
                runNames = [runNames {sprintf('run %d', kk)}];
            end
            plot(1:nStim, squeeze(mean(ts(el,:,:),3)), 'k', 'LineWidth', 2); axis tight;
            runNames = [runNames {'average'}];
            title(plotTitle);
            if el == 1; xlabel('PRF stimulus (#)'); ylabel('Broadband signal change'); end %legend(runNames); 
            %set(gcf, 'Position', [150 100 1500 1250]);
            set(gcf, 'Position', get(0,'screensize'));
            set(gca, 'FontSize', 14);
        end
        saveas(gcf, fullfile(plotSaveDir,'data', figureName), 'png'); close;        
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
	channels = fulldata{ii}.channels;
    
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
            save(saveName, 'channels','data','results','tr','opt','subject');  
        end
    
        % Make plots of the estimated PRFs and PRF fits

        if doPlots

            % Timeseries + fits
            figureName = sprintf('%s_prfmodelfits', subject);
            figure('Name', figureName); plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
            
            ecog_plotPRFtsfits(data,results,tr, opt,channels,subject)
            
            saveas(gcf, fullfile(plotSaveDir,'modelfits', figureName), 'png'); close;

            % PRFs
            coloropt = 0;
            figureName = sprintf('%s_prfs', subject);
            figure('Name', figureName); plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
            
            ecog_plotPRFs(channels, results, stimulus, coloropt)  
            
            saveas(gcf, fullfile(plotSaveDir,'modelfits', figureName), 'png'); close;
            
            % Add some fun summary plots comparing e.g. benson and
            % analyzePRF?
        end
    end
end