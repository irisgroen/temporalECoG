function [data] = tde_selectData(data, stimNames, epochOpts, elecOpts, baselineTime, savePlots, plotSaveDir)

% Description
%
% function [data] = tde_selectData(data, stimNames, epochOpts, elecOpts, baselineTime)
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

% <stimNames>
if ~exist('stimNames', 'var') || isempty(stimNames)
    stimNames = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5',...
                 'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
                 'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};
end

% <epochOpts>
if ~exist('epochOpts', 'var') || isempty(epochOpts)
    epochOpts = struct(); % a struct specifying criteria for inclusion of epochs
    epochOpts.outlier_thresh= 5; % x-fold magnitude above which epoch will be labeled as outlier
    epochOpts.maxnooutlier  = 50; % number of outlier epochs after which entire channel will be labeled as bad
    epochOpts.write         = 1; % boolean
end

% <elecOpts>
if ~exist('elecOpts', 'var') || isempty(elecOpts)
    elecOpts = struct(); % a struct specifying criteria for inclusion of electrodes
    elecOpts.max_thresh     = 1; % minimum required maximal response in % signal change
    elecOpts.mean_thresh    = 0; % minimum required mean response during stim_on period in % signal change
    elecOpts.stim_on        = [0 0.2]; % period across which to take mean
    elecOpts.exclude_bad    = 1; % boolean
    elecOpts.exclude_depth  = 1; % boolean
end

% <baselineTime>
if ~exist('baselineTime', 'var') || isempty(baselineTime)
    baselineTime = [-0.2 0];
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

%%

% Loop over subjects
for ii = 1:length(data)
    
    subject  = data{ii}.name;
    epochs   = data{ii}.epochs;
    channels = data{ii}.channels;
    events   = data{ii}.events;
    t        = data{ii}.t;
    
    %% STEP 1 CHECK FOR OUTLIERS EPOCHS/CHANNELS
%     fprintf('[%s] Checking for bad epochs ...\n',mfilename);
%     % remove epochs that have 20x the average
%     % if more than X epochs for a channel, remove entire channel
%     % output a description of how many trials were removed (write to
%     % file?)
%     % update the events file? (means we're removing epoch from all chans)
%     
%     % compute sum over time for each epoch
%     summed_epochs = squeeze(sum(epochs,2));
%     %figure;histogram(summed_epochs(1,:),100)
% 
%     chan_mean = mean(summed_epochs,2);
%     for jj = 1:size(chan_mean,1)
%         %figure;histogram(summed_epochs(jj,:),100)
%         outlier_idx1 = find(summed_epochs(jj,:) > chan_mean(jj) + (epochOpts.outlier_thresh-1) * chan_mean(jj));
%         outlier_idx2 = find(summed_epochs(jj,:) < chan_mean(jj) - epochOpts.outlier_thresh * chan_mean(jj));
%         outlier_idx = [outlier_idx1 outlier_idx2];
%         %keep_idx = setdiff(1:size(epochs,3),outlier_idx);
%         disp(length(outlier_idx));
%         epochs(jj,:,outlier_idx) = nan; 
%     end
       

  %% STEP 2 convert to percent signal change 
    fprintf('[%s] Converting epochs to percent signal change...\n',mfilename);

    % Provide run index to perform separately for each run and session
    [~,~,ses_idx]= unique(events.session_name);
    [~,~,run_idx] = unique(events.run_name);
    idx = (ses_idx*100)+run_idx;

    [epochs] = ecog_normalizeEpochs(epochs, t, baselineTime, 'percentsignalchange', idx);
    channels.units = repmat({'%change'}, [height(channels),1]);
    
 %% STEP 3 select electrodes
        
    % Restrict selection to relevant stimuli only
    stimsForSelection = contains(events.trial_name, stimNames);
  
    mean_resp = mean(epochs(:,:,stimsForSelection),3);
    llim = (mean_resp - (std(epochs(:,:,stimsForSelection),0,3)));
    ulim = (mean_resp + (std(epochs(:,:,stimsForSelection),0,3)));
    mean_resp_ci = cat(3, llim, ulim);
    nEl = size(mean_resp,1); 

	% plot
    if savePlots            
        figureName = sprintf('viselec_%s_all', subject);
        figure('Name', figureName); nSubPlot = ceil(sqrt(nEl)); 
        for el = 1:nEl
            subplot(nSubPlot,nSubPlot,el); hold on
            plotTitle = sprintf('%s W:%s B:%s ', channels.name{el}, channels.wangarea{el}, channels.bensonarea{el});        
            ecog_plotSingleTimeCourse(t, mean_resp(el,:), squeeze(mean_resp_ci(el,:,:)), [], plotTitle);
            if el == 1; xlabel('Time (s)'); ylabel('Broadband signal change');end
            set(gcf, 'Position', [150 100 1500 1250]);
        end
        saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
    end
    
    stimOnIdx = t > elecOpts.stim_on(1) & t <= elecOpts.stim_on(2);  
    % FIRST SELECTION CRITERION:
    select_idx1 = max(mean_resp(:,stimOnIdx),[],2) > elecOpts.max_thresh;

    % SECOND SELECTION CRITERION:
    select_idx2 = mean(mean_resp(:,stimOnIdx),2) > elecOpts.mean_thresh;
    %select_idx2 = zeros(size(select_idx1));
    %for el = 1:nEl
    %    if mean(mean_resp(el, elecOpts.stim_on)) > elecOpts.max_thresh; select_idx2(el) = 1; end
    %end

    % EXCLUDE CHANNELS LABELED AS BAD
	select_idx3 = ones(size(select_idx1));
    if elecOpts.exclude_bad, select_idx3 = contains(channels.status, 'good');end  

    % EXCLUDE DEPTH ELECTRODES
    select_idx4 = ones(size(select_idx1));
    if elecOpts.exclude_depth, select_idx4 = contains(lower(channels.type), 'ecog');end  

    % Combine criteria
    select_idx = select_idx1 + select_idx2 + select_idx3 + select_idx4;
    select_idx = (select_idx == 4);

    % Plot selection
    if savePlots           
        figureName = sprintf('viselec_%s_selected', subject);
        figure('Name', figureName); nSubPlot = ceil(sqrt(nEl)); 
        for el = 1:nEl
            if select_idx(el)
                subplot(nSubPlot,nSubPlot,el); hold on
                plotTitle = sprintf('%s W:%s B:%s ', channels.name{el}, channels.wangarea{el}, channels.bensonarea{el});
                ecog_plotSingleTimeCourse(t, mean_resp(el,:), squeeze(mean_resp_ci(el,:,:)), [], plotTitle)    
                set(gcf, 'Position', [150 100 1500 1250]);
            end
        end
        saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
    end
    
    % Select the channels, update channels table
    epochs = epochs(select_idx,:,:);
    channels = channels(select_idx,:);
    
 %% STEP 4 average across conditions?
 
 %% STEP 5 plot selection? (single trials, averages?)
 
        % Checks (Temporary; TO write simpler function e.g. bidsEcogPlotEpochs)
%         trials = [];
%         trials.broadband = epochs;
%         trials.events = events;
%         trials.bb_bands = [50 200];
%         trials.time = t;
%         trials.channels = channels;
%         trials.viselec = visualelectrodes;
%         whichElectrodes = channels.name;
%         specs.baselineType = 'selectedtrials';
%         specs.dataTypes = {'broadband'};
%         trialType = {'ONEPULSE-1', 'ONEPULSE-2', 'ONEPULSE-3', 'ONEPULSE-4', 'ONEPULSE-5', 'ONEPULSE-6'};
%         ecog_plotTimecourses(trials, whichElectrodes, trialType, specs);
%         trialType = {'TWOPULSE-1', 'TWOPULSE-2', 'TWOPULSE-3', 'TWOPULSE-4', 'TWOPULSE-5', 'TWOPULSE-6'};
%         ecog_plotTimecourses(trials, whichElectrodes, trialType, specs);
%         trialType = {'CRF-1', 'CRF-2', 'CRF-3', 'CRF-4', 'CRF-5'};
%         ecog_plotTimecourses(trials, whichElectrodes, trialType, specs);

    %% STEP 6 update channel and event tables, save?
    
	data{ii}.epochs     = epochs;
    data{ii}.channels   = channels;
    data{ii}.events     = events;
    
end