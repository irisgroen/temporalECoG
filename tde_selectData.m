function [allData, allChannels, stimNames, t] = tde_selectData(data, savePlots, plotSaveDir, stimNames, epochOpts, elecOpts, baselineTime)

% Description
%
% function [allData, allChannels, stimNames, t] = tde_selectData(data_in, ...
%           [savePlots], [plotSaveDir], [stimNames], [epochOpts], [elecOpts], [baselineTime])
% 
% Removes bad epochs and channels with many bad epochs (epochOpts)
% Removes channels that do not match inclusion criteria (elecOpts)
% Converts to percent signal change (using baselineTime)
% Outputs reduced version of data 
% Makes plots

% <data>
if ~exist('data', 'var') || isempty(data)
	error('Please provide the data struct outputted by tde_getData,m as input');
end 

% <makePlots>
if ~exist('savePlots', 'var') || isempty(savePlots)
    savePlots = 1;
end

% <plotSaveDir>
if ~exist('plotSaveDir', 'var') || isempty(plotSaveDir)
    plotSaveDir = '/Volumes/server/Projects/BAIR/Papers/TemporalDynamicsECoG/figures/electrodeselection';
    if ~exist(plotSaveDir, 'dir'); mkdir(plotSaveDir); end
end

% <stimNames>
if ~exist('stimNames', 'var') || isempty(stimNames)
    stimNames = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5',...
                 'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
                 'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};
end

% <epochOpts>
if ~exist('epochOpts', 'var') || isempty(epochOpts)
    epochOpts = struct(); % a struct specifying criteria for inclusion of epochs
    epochOpts.outlier_thresh= 20; % x-fold magnitude above which epoch will be labeled as outlier
    %epochOpts.maxnooutlier  = 50; % number of outlier epochs after which entire channel will be labeled as bad
    %epochOpts.write         = 1; % boolean
end

% <elecOpts>
if ~exist('elecOpts', 'var') || isempty(elecOpts)
    elecOpts = struct(); % a struct specifying criteria for inclusion of electrodes
    elecOpts.max_thresh     = 1; % minimum required maximal response in % signal change
    elecOpts.mean_thresh    = 0; % minimum required mean response during stim_on period in % signal change
    elecOpts.stim_on        = [0 0.5]; % time period across which to compute mean response
    elecOpts.exclude_bad    = 1; % boolean
    elecOpts.exclude_depth  = 0; % boolean
end

% <baselineTime>
if ~exist('baselineTime', 'var') || isempty(baselineTime)
    baselineTime = [-0.2 0];
end

%% Loop over subjects

% Initialize
allData    = [];
allChannels = [];
    
for ii = 1:length(data)
    
    subject  = data{ii}.subject;
    epochs   = data{ii}.epochs;
    channels = data{ii}.channels;
    events   = data{ii}.events;
    t        = data{ii}.t;
    
    fprintf('[%s] Selecting data for subject %s \n',mfilename, subject);

    %% STEP 1 CHECK FOR OUTLIERS EPOCHS/CHANNELS
    %fprintf('[%s] Removing bad epochs...\n',mfilename);
    
    % Remove epochs whose summed amplitude is outlier thresh x more or less
    % than the average across all channels
    
    % Compute sum over time for each epoch
    summed_epochs = squeeze(sum(epochs,1));
    newepochs = epochs;
    for jj = 1:height(channels)
        outlier_idx1 = summed_epochs(:,jj) > epochOpts.outlier_thresh * median(summed_epochs(:,jj));
        outlier_idx2 = summed_epochs(:,jj) < summed_epochs(:,jj) - (epochOpts.outlier_thresh * summed_epochs(:,jj));
        outlier_idx = logical(outlier_idx1+outlier_idx2);
        newepochs(:,outlier_idx, jj) = nan;
        
        outliers = find(outlier_idx);
        if ~isempty(outliers)
            figureName = sprintf('outlierepochs_sub-%s_chan-%s', subject, channels.name{jj});
            figure('Name', figureName); hold on;
            nOutliers = length(outliers);
            dim1 = round((nOutliers+1)/2);
            dim2 = round((nOutliers+1)/dim1);
            subplot(dim1,dim2,1); hold on; 
            histogram(summed_epochs(:,jj),100); title(channels.name{jj});
            for kk = 1:nOutliers
                subplot(dim1,dim2,kk+1); %plot(t, epochs(:,outliers(kk), jj), 'LineWidth');
                ecog_plotSingleTimeCourse(t, epochs(:,outliers(kk),jj), [], [], sprintf('epoch %d', outliers(kk)));    
            end
            set(gcf, 'Position', [150 100 1500 1250]);
            saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
        end
    end
        
    % TO DO include plots of single trials pre and post removal
    % if more than X epochs for a channel, remove entire channel
    % output a description of how many trials were removed (write to
    % file?), also update events file?    
    
  %% STEP 2 convert to percent signal change 
    %fprintf('[%s] Converting epochs to percent signal change...\n',mfilename);

    % Provide run index to perform separately for each run and session
    [~,~,task_idx]= unique(events.task_name);
    [~,~,ses_idx]= unique(events.session_name);
    [~,~,run_idx] = unique(events.run_name);
    [~,~,idx] = unique([task_idx ses_idx run_idx], 'rows');
    idx = [];
    [epochs] = ecog_normalizeEpochs(epochs, t, baselineTime, 'percentsignalchange', idx);
    channels.units = repmat({'%change'}, [height(channels),1]);
    
  %% STEP 3 select electrodes   
    %fprintf('[%s] Selecting electrodes...\n',mfilename);

    % Restrict selection to relevant stimuli only
    stimsForSelection = contains(events.trial_name, stimNames);
  
    mean_resp = mean(epochs(:,stimsForSelection,:),2, 'omitnan');
    llim = (mean_resp - (std(epochs(:,stimsForSelection,:),0,2,'omitnan')));
    ulim = (mean_resp + (std(epochs(:,stimsForSelection,:),0,2,'omitnan')));
    mean_resp_sd = cat(2, llim, ulim);

	% plot
    if savePlots            
        nEl = size(mean_resp,3); 
        figureName = sprintf('viselec_%s_all', subject);
        figure('Name', figureName); plotDim1 = round(sqrt(nEl)); plotDim2 = ceil((nEl)/plotDim1);
        for el = 1:nEl
            subplot(plotDim1,plotDim2,el); hold on
            plotTitle = sprintf('%s W:%s B:%s ', channels.name{el}, channels.wangarea{el}, channels.bensonarea{el});        
            ecog_plotSingleTimeCourse(t, mean_resp(:,:,el), squeeze(mean_resp_sd(:,:,el)), [], plotTitle);
            %if el == 1; xlabel('Time (s)'); ylabel('Broadband signal change');end
            set(gcf, 'Position', [150 100 1500 1250]);
        end
        saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
    end
    
    select_idx = ones(height(channels),4);
    
    % Exclude channels based on max and mean:
    stim_on_idx = t > elecOpts.stim_on(1) & t <= elecOpts.stim_on(2);  
    select_idx(:,1) = max(mean_resp(stim_on_idx,:),[],1) > elecOpts.max_thresh;
    select_idx(:,2) = mean(mean_resp(stim_on_idx,:),1) > elecOpts.mean_thresh;

    % EXCLUDE CHANNELS LABELED AS BAD
    if elecOpts.exclude_bad, select_idx(:,3) = contains(channels.status, 'good');end  

    % EXCLUDE DEPTH ELECTRODES
    if elecOpts.exclude_depth, select_idx(:,3) = contains(lower(channels.type), 'ecog');end  

    % Combine criteria
    select_idx = sum(select_idx,2) == 4;

    % Plot selection
    if savePlots           
        figureName = sprintf('viselec_%s_selected', subject);
        figure('Name', figureName); 
        for el = 1:nEl
            if select_idx(el)
                subplot(plotDim1,plotDim2,el); hold on
                plotTitle = sprintf('%s W:%s B:%s ', channels.name{el}, channels.wangarea{el}, channels.bensonarea{el});
                ecog_plotSingleTimeCourse(t, mean_resp(:,el), squeeze(mean_resp_sd(:,:,el)), [], plotTitle)    
                set(gcf, 'Position', [150 100 1500 1250]);
            end
        end
        saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
    end
    
    % Select the channels
    epochs = epochs(:,:,select_idx);
    % Update channels table
    channels = channels(select_idx,:);
    
    %% STEP 4 average across trials, concatenate subjects
    
    % Average across trials within stimulus condition
    epochs_averaged = nan(size(epochs,1), length(stimNames), size(epochs,3));
    for jj = 1:length(stimNames)
        trial_idx = contains(events.trial_name, stimNames{jj});
        epochs_averaged(:,jj,:) = nanmean(epochs(:,trial_idx,:),2);
    end
    % Also compute se across trials?
    
    % Concatenate the data across subjects
    allData = cat(3, allData, epochs_averaged);    
    
    % Remove a number of columns from channel table for readability, and
    % concatenate across subjects
    channels = removevars(channels, {'low_cutoff', 'high_cutoff', 'reference', 'group', 'sampling_frequency', 'bb_method', 'bb_bandwidth', 'status'});
    allChannels = [allChannels; channels];
    
    % TO DO: plot with final time courses per condition for each sub
    
end
fprintf('[%s] Done! \n',mfilename);
end