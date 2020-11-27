function [data] = tde_computePRFtimecourses(data, time_win, normalizeEpochs, doPlots, plotSaveDir)

% Computes prf timecourses for data to be fitted with analyzePRF
% [data2fit] = tde_computePRFs(recomputeData, doPlots, saveDir, resultsStr) 
%
% <data> data struct computed by tde_getData.m
% <timeWin> time window over which to average the broadband timecourse to
%   get the PRF activity estimate for each bar position
%   default: [0.05 0.55];
% <doPlots> flag indicating whether to save out plots of data and
%   prfs to fullfile(analysisRootPath, 'figures', 'prfs')
%   default 'false

% <time_win>
if ~exist('time_win','var') || isempty(time_win)
    time_win  = [0.05 0.55];
end

% <normalize>
if ~exist('normalizeEpochs','var') || isempty(normalizeEpochs)
    normalizeEpochs = false;
end

% <doPlots>
if ~exist('doPlots','var') || isempty(doPlots)
    doPlots = false; % boolean
end

% <plotSaveDir>
if ~exist('plotSaveDir','var') || isempty(plotSaveDir)
    plotSaveDir = fullfile(analysisRootPath, 'figures', 'prfs');
end
if ~exist(plotSaveDir, 'dir'); mkdir(fullfile(plotSaveDir));end

% Compute PRF timecourses for each subject
nSubjects = length(data);

for ii = 1:nSubjects
    
    subject  = data{ii}.subject;
    t        = data{ii}.t;
    epochs   = data{ii}.epochs_b;
    
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
    
    % Determine number of runs and stimuli for this subject
    [~, nTrials, nChans] = size(epochs);
    nStim = length(stimInx);    
    nRuns = nTrials/nStim; 
       
    % Normalize epochs, within run?
    if normalizeEpochs
        run_indices = []; 
        for jj = 1:nRuns
            run_indices = [run_indices; ones(nStim,1)*jj];
        end
        [epochs] = ecog_normalizeEpochs(epochs, t, [], [], run_indices);
    end
    
    % Compute average broadband response in time window
    trials = squeeze(mean(epochs(t>time_win(1) & t<time_win(2),:,:),1)); 

    % Transpose to have channels in first dimension
    trials = trials';
    
    % Reshape to separate individual runs    
    ts = reshape(trials,[nChans nStim nRuns]);
    
    data{ii}.ts       = ts;
    data{ii}.stim_inx = stimInx;
    
    % Make plots of the trials and of the PRF timecourses
    
    if doPlots
        
        channels = data{ii}.channels;

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
            %if el == 1; xlabel('Time (s)'); ylabel('Broadband signal change');end
            set(gca, 'XTickLabel', []);
        end
        %set(gcf, 'Position', [150 100 1500 1250]);
        set(gcf, 'Position', get(0,'screensize'));
        set(findall(gcf,'-property','FontSize'),'FontSize',14)
        saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
        
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
            set(gca, 'XTickLabel', []);
            %if el == 1; xlabel('PRF stimulus (#)'); ylabel('Broadband signal change'); end %legend(runNames); 
        end
        %set(gcf, 'Position', [150 100 1500 1250]);
        set(gcf, 'Position', get(0,'screensize'));
        set(findall(gcf,'-property','FontSize'),'FontSize',14)
        saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;        
    end   
end
   