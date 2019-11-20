function [data, channels, stimNames, t, srate] = tde_selectData(inputData, stimNames, opts)

% Description
%
% function [data, channels, stimNames, t] = tde_selectData(inputData, [stimNames], [opts])
% 
% Outputs reduced version of data after following steps:
%
% Removes bad epochs and channels with many bad epochs (epochOpts)
% Removes channels that do not match inclusion criteria (elecOpts)
% Converts to percent signal change (using baselineTime)
% Averages across trials (make optional?)
% Normalizes by max (optional)
% Averages across visual areas (optional)
% Makes plots (optional)

%% Check inputs

% <data>
if ~exist('inputData', 'var') || isempty(inputData)
	error('Please provide the data struct outputted by tde_getData.m as input');
end 

% <stimNames>
if ~exist('stimNames', 'var') || isempty(stimNames)
    stimNames = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5',...
                 'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
                 'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};
end

% <opts>
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end
if ~isfield(opts,'stim_on') || isempty(opts.stim_on)
    opts.stim_on             = [0 0.5]; % time period across which stimulus is presented
end
if ~isfield(opts,'baseline_time') || isempty(opts.baseline_time)
    opts.baseline_time       = [-0.2 0]; % time period across which to compute normalization baseline
end
if ~isfield(opts,'epoch_outlier_thresh') || isempty(opts.epoch_outlier_thresh)
    opts.epoch_outlier_thresh= 20; % x-fold max magnitude above which epoch will be labeled as outlier
end
if ~isfield(opts,'elec_max_thresh') || isempty(opts.elec_max_thresh)
    opts.elec_max_thresh     = 1; % minimum required maximal response in % signal change for electrode inclusion
end
if ~isfield(opts,'elec_mean_thresh') || isempty(opts.elec_mean_thresh)
    opts.elec_mean_thresh    = 0; % minimum required mean response during stim_on period in % signal change
end
if ~isfield(opts,'elec_exclude_depth') || isempty(opts.elec_exclude_depth)
    opts.elec_exclude_depth  = false; % boolean
end
if ~isfield(opts,'average_trials') || isempty(opts.average_trials)
    opts.average_trials      = true;  % boolean
end
if ~isfield(opts,'normalize_data') || isempty(opts.normalize_data)
    opts.normalize_data      = true;  % boolean
end
if ~isfield(opts,'average_elecs') || isempty(opts.average_elecs)
    opts.average_elecs       = false; % boolean
end
if ~isfield(opts,'doplots') || isempty(opts.doplots)
    opts.doplots             = true; % boolean
end
if ~isfield(opts,'plotsavedir') || isempty(opts.plotsavedir)
    opts.plotsavedir         = '/Volumes/server/Projects/BAIR/Papers/TemporalDynamicsECoG/figures/electrodeselection';
end

savePlots   = opts.doplots;
plotSaveDir = opts.plotsavedir;
if ~exist(plotSaveDir, 'dir'); mkdir(plotSaveDir); end

%% Loop over subjects

% Initialize
data     = [];
channels = [];
    
for ii = 1:length(inputData)
    
    subject     = inputData{ii}.subject;
    epochs      = inputData{ii}.epochs;
    channels_in = inputData{ii}.channels;
    events      = inputData{ii}.events;
    t           = inputData{ii}.t;

    fprintf('[%s] Selecting data for subject %s \n',mfilename, subject);
          
    % Restrict selection to relevant stimuli only
    stimsForSelection = contains(events.trial_name, stimNames);
    epochs = epochs(:, stimsForSelection, :);
    events = events(stimsForSelection, :);
    
    % Compute time indices
    stim_on_idx = t > opts.stim_on(1) & t <= opts.stim_on(2);  

    %% STEP 1 CHECK FOR OUTLIERS EPOCHS/CHANNELS
    
    %fprintf('[%s] Removing bad epochs...\n',mfilename);
    
    % Remove epochs whose max response amplitude is outlier thresh x more
    % or less than the average response in that channel
    
    % Compute max over time for each epoch
    max_epochs_alltmpts = squeeze(max(epochs,[],1));
    max_epochs_stim_on = squeeze(max(epochs(stim_on_idx,:,:),[],1));
    for jj = 1:height(channels_in)
        outlier_thresh = opts.epoch_outlier_thresh * median(max_epochs_stim_on(:,jj));
        outlier_idx = max_epochs_alltmpts(:,jj) > outlier_thresh;             
        outliers = find(outlier_idx);
        if savePlots            
            if ~isempty(outliers)
                figureName = sprintf('outlierepochs_sub-%s_chan-%s', subject, channels_in.name{jj});
                figure('Name', figureName); hold on;
                nOutliers = length(outliers);
                dim1 = round((nOutliers+1)/2);
                dim2 = round((nOutliers+1)/dim1);
                subplot(dim2,dim1,1); hold on; title(channels_in.name{jj}); 
                histogram(max_epochs_alltmpts(:,jj),100); line([outlier_thresh outlier_thresh], get(gca, 'YLim'), 'Color', 'r','LineStyle', ':', 'LineWidth', 2);
                set(gca, 'fontsize', 14); xlabel('max broadband'); ylabel('number of epochs');
                for kk = 1:nOutliers
                    subplot(dim2,dim1,kk+1); 
                    ecog_plotSingleTimeCourse(t, epochs(:,outliers(kk),jj), [], [], sprintf('epoch %d %s', outliers(kk), events.trial_name{outliers(kk)}));    
                end
                set(gcf, 'Position', [150 100 300*dim1 300*dim2]);
                saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
            end
        end
        epochs(:,outlier_idx,jj) = nan;
    end
    
  %% STEP 2 convert to percent signal change 
    %fprintf('[%s] Converting epochs to percent signal change...\n',mfilename);

    % Provide run index to perform separately for each run and session
%     [~,~,task_idx]= unique(events.task_name);
%     [~,~,ses_idx]= unique(events.session_name);
%     [~,~,run_idx] = unique(events.run_name);
%     [~,~,idx] = unique([task_idx ses_idx run_idx], 'rows');
    idx = [];
    [epochs] = ecog_normalizeEpochs(epochs, t, opts.baseline_time, 'percentsignalchange', idx);
    channels_in.units = repmat({'%change'}, [height(channels_in),1]);
    
  %% STEP 3 select electrodes   
    %fprintf('[%s] Selecting electrodes...\n',mfilename);

    mean_resp = mean(epochs,2, 'omitnan');
    llim = (mean_resp - (std(epochs,0,2,'omitnan')));
    ulim = (mean_resp + (std(epochs,0,2,'omitnan')));
    mean_resp_sd = cat(2, llim, ulim);

	% plot
    if savePlots            
        nEl = size(mean_resp,3); 
        figureName = sprintf('viselec_%s_all', subject);
        figure('Name', figureName); plotDim1 = round(sqrt(nEl)); plotDim2 = ceil((nEl)/plotDim1);
        for el = 1:nEl
            subplot(plotDim1,plotDim2,el); hold on
            plotTitle = sprintf('%s W:%s B:%s ', channels_in.name{el}, channels_in.wangarea{el}, channels_in.bensonarea{el});        
            ecog_plotSingleTimeCourse(t, mean_resp(:,:,el), squeeze(mean_resp_sd(:,:,el)), [], plotTitle);
            %if el == 1; xlabel('Time (s)'); ylabel('Broadband signal change');end
            set(gcf, 'Position', [150 100 1500 1250]);
        end
        saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
    end
    
    select_idx = ones(height(channels_in),3);
    
    % Exclude channels based on max and mean:
    stim_on_idx = t > opts.stim_on(1) & t <= opts.stim_on(2);  
    select_idx(:,1) = max(mean_resp(stim_on_idx,:),[],1) > opts.elec_max_thresh;
    select_idx(:,2) = mean(mean_resp(stim_on_idx,:),1) > opts.elec_mean_thresh;

    % Exclude depth electrodes
    if opts.elec_exclude_depth, select_idx(:,3) = contains(lower(channels_in.type), 'ecog');end  

    % Combine criteria
    select_idx = sum(select_idx,2) == 3;

    % Plot selection
    if savePlots           
        figureName = sprintf('viselec_%s_selected', subject);
        figure('Name', figureName); 
        for el = 1:nEl
            if select_idx(el)
                subplot(plotDim1,plotDim2,el); hold on
                plotTitle = sprintf('%s W:%s B:%s ', channels_in.name{el}, channels_in.wangarea{el}, channels_in.bensonarea{el});
                ecog_plotSingleTimeCourse(t, mean_resp(:,el), squeeze(mean_resp_sd(:,:,el)), [], plotTitle)    
                set(gcf, 'Position', [150 100 1500 1250]);
            end
        end
        saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
    end
    
    % Select the channels
    epochs = epochs(:,:,select_idx);
    % Update channels table
    channels_in = channels_in(select_idx,:);
    
    %% STEP 4 average across trials, concatenate subjects
    
    % Average trials
    if opts.average_trials
        [epochs] = average_trials(epochs, events, stimNames);
    end
    
    % Concatenate the data across subjects
    data = cat(3, data, epochs);    
    
    % Remove a number of columns from channel table for readability, and
    % concatenate across subjects
    channels_in = removevars(channels_in, {'low_cutoff', 'high_cutoff', 'reference', 'group', 'bb_method', 'bb_bandwidth'});
    channels = [channels; channels_in];
    
end

%% Run checks and finalize processing %%

% Now that we have the selected trials and electrodes from all subjects, 
% run a few checks, and make decisions on how to further process the data.

% Check that all channels have some sample frequencies
assert(length(unique(channels.sampling_frequency))==1);
srate = channels.sampling_frequency(1);

% Sort electrodes on visual area (rather than subjectID)?
if opts.sort_channels
    [channels] = sort_channels(channels);
end

% Add an index column to channels
index = [1:height(channels)]';
channels = addvars(channels, index, 'Before', 'name');

% Scale each electrode to its max?
if opts.normalize_data       
    [data] = normalize_data(data);
end

% Average elecs within area?
if opts.average_elecs
    [data, channels] = average_elecs(data, channels);
end

% Plot final selected data
%if opts.doplots
    figure, sz = ceil(sqrt(size(data,3)));
    for ii = 1:size(data,3)
        subplot(sz,sz,ii); plot(data(:,:,ii), 'LineWidth', 3); 
        title(sprintf('%d %s:%s B:%s W:%s', ...
            channels.index(ii), channels.subject_name{ii}, channels.name{ii}, channels.bensonarea{ii}, channels.wangarea{ii})); 
    end
    set(gcf, 'Position', [400 200 2000 1200]);
%end

fprintf('[%s] Done! \n',mfilename);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Trial averaging
function [epochs_averaged] = average_trials(epochs, events, stimNames)   
    
    % Average across trials within stimulus condition
    epochs_averaged = nan(size(epochs,1), length(stimNames), size(epochs,3));
    for ii = 1:length(stimNames)
        trial_idx = contains(events.trial_name, stimNames{ii});
        epochs_averaged(:,ii,:) = nanmean(epochs(:,trial_idx,:),2);
    end
    % Also compute se across trials?
end

% Channel sorting
function [channels] = sort_channels(channels)   
    % sort on benson area
    [~,I] = sort(channels.bensonarea);
    channels = channels(I,:);
end

% Data normalization
function [data] = normalize_data(data)
    normdata = data;
    % normalize each channel separately
    for ii = 1:size(data,3)
        tmp = data(:, :, ii);
        maxRsp(ii) = max(tmp(:));
        normdata(:, :, ii) = data(:, :, ii)./maxRsp(ii);
    end
    data = normdata;
end

% Data averaging
function [data, channels] = average_elecs(data, channels)
    
    avdata = nan(size(data,1), size(data,2), 4);

    INX = [];
    INX{1} = contains(channels.wangarea, 'V1') | contains(channels.bensonarea, 'V1');
    INX{2} = contains(channels.wangarea, 'V2') | contains(channels.bensonarea, 'V2');
    INX{3} = contains(channels.wangarea, 'V3') & ~contains(channels.bensonarea, {'V3a', 'V3b'}) | contains(channels.bensonarea, 'V3') & ~contains(channels.bensonarea, {'V3a', 'V3b'});
    INX{4} = contains(channels.wangarea, {'V3a', 'V3b', 'IPS', 'hV4', 'LO', 'TO'}) | contains(channels.bensonarea, {'V3a', 'V3b', 'hV4', 'LO', 'TO'}) ;

    for ii = 1:length(INX)
        mdata = mean(data(:,:,INX{ii}),3);
        
        if opts.normalize_data     
            % another normalization step from Jings code, necessary?
            maxmdata = max(mdata(:));
            mdata    = mdata ./maxmdata;
        end
        
        avdata(:,:,ii) = mdata;
    end
    data = avdata;
    channels = []; 
    
    % TO DO prepare a new channel table with {'V1', 'V2', 'V3', 'Vhigher'};
end
end