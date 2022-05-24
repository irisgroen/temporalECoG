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
%
% 2020 Iris Groen

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
        case 'p01'
            % in this subject, trigger 41 was not sent because it was the
            % same as used for the blank (stimulus coding error).
            stimInx = setdiff(1:224, 41);
        case 'p04'
            % In this subject, a different version of the PRF experiment
            % was run, with different apertures. Since this subject only
            % contributes a few channels, we'll skip their PRF analysis.
            continue
        case 'p05'
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
        f_ind = checkForHDgrid(channels);
        
        nEpochs  = size(epochs,2);
        
        % Plot individual trials
        fprintf('[%s] Plotting prf trials for subject %s \n',mfilename, subject);
        for f = 1:length(f_ind)
            if length(f_ind) > 1
                figureName = sprintf('%s_prftrials_%d', subject, f);
            else
                figureName = sprintf('%s_prftrials', subject);
            end
            nChans = length(f_ind{f});
            figure('Name', figureName); plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
            for el = 1:nChans
                subplot(plotDim1,plotDim2,el); hold on
                el_ind = f_ind{f}(el);
                for kk = 1:nEpochs
                    plot(t, epochs(:,kk,el_ind)); 
                end
                plot(t, squeeze(mean(epochs(:,kk,el_ind),2)), 'k', 'LineWidth', 3); axis tight;

                yLim = get(gca, 'YLim');
                line([0 0], yLim,'LineStyle', ':', 'Color', 'k');
                line([t(1) t(end)], [0 0],'LineStyle', ':', 'Color', 'k');
                plotTitle = sprintf('%s %s %s ', channels.name{el_ind}, channels.benson14_varea{el_ind}, channels.wang15_mplbl{el_ind});        
                title(plotTitle);
                axis tight
                set(gca, 'XTickLabel', []);
            end
            set(gcf, 'Position', get(0,'screensize'));
            set(findall(gcf,'-property','FontSize'),'FontSize',14)
            saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
        end
        
        % Plot PRF timecourses for each run + average
        fprintf('[%s] Plotting prf timecourses data for subject %s \n',mfilename, subject);
        for f = 1:length(f_ind)
            if length(f_ind) > 1
                figureName = sprintf('%s_prftimecourses_%d', subject, f);
            else
                figureName = sprintf('%s_prftimecourses', subject);
            end
            nChans = length(f_ind{f});

            figure('Name', figureName); plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
            for el = 1:nChans
                subplot(plotDim1,plotDim2,el); hold on
                el_ind = f_ind{f}(el);
                plotTitle = sprintf('%s %s %s ', channels.name{el_ind}, channels.benson14_varea{el_ind}, channels.wang15_mplbl{el_ind});        
                runNames = [];
                for kk = 1:nRuns
                    plot(1:nStim, squeeze(ts(el_ind,:,kk)), 'LineWidth', 1);
                    runNames = [runNames {sprintf('run %d', kk)}];
                end
                plot(1:nStim, squeeze(mean(ts(el_ind,:,:),3)), 'k', 'LineWidth', 2); axis tight;
                runNames = [runNames {'average'}];
                title(plotTitle);
                set(gca, 'XTickLabel', []);
                set(gca, 'YTickLabel', []);

            end
            set(gcf, 'Position', get(0,'screensize'));
            set(findall(gcf,'-property','FontSize'),'FontSize',14)
            saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;  
        end
    end   
end
   