function [epochs, channels, stimnames, t, srate, opts, epochs_se, epochs_indiv] = tde_selectData(data, opts)

% Description
%
% [epochs, channels, stimnames, t, srate, ...
%   opts, epochs_se, epochs_indiv] = tde_selectData(data, [opts])
% 
% Outputs reduced version of data after following steps:
%
% Removes bad epochs and channels with many bad epochs (opts.epoch_)
% Removes channels that do not match inclusion criteria (opts.elec_)
% Converts to percent signal change (using baselineTime)
% Averages across trials (make optional?)
% Normalizes by max (optional)
% Averages across visual areas (optional)
% Makes plots (optional)

%% Check inputs

% <data>
if ~exist('data', 'var') || isempty(data)
	error('Please provide the data struct outputted by tde_getData.m as input');
end 

% <opts>
if ~exist('opts','var') || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 'stimnames') || isempty(opts.stimnames)
    opts.stimnames = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5',...
                 'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
                 'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};
end
if ~isfield(opts, 'areanames') 
    opts.areanames = [];
end
if ~isfield(opts,'stim_on') || isempty(opts.stim_on)
    opts.stim_on             = [0 1]; % time period across which stimulus is presented
end
if ~isfield(opts,'baseline_time') || isempty(opts.baseline_time)
    opts.baseline_time       = [-0.2 0]; % time period across which to compute normalization baseline
end
if ~isfield(opts,'epoch_jump_thresh') || isempty(opts.epoch_jump_thresh)
    opts.epoch_jump_thresh = 500; % max jump in voltage allowed within stim_on period
end
if ~isfield(opts,'epoch_outlier_thresh') || isempty(opts.epoch_outlier_thresh)
    opts.epoch_outlier_thresh = 20; % percentile for max absolute voltage within stim_on period at which epoch will be labeled as outlier
end
if ~isfield(opts,'elec_selection_method') || isempty(opts.elec_selection_method)
    opts.elec_selection_method = 'splithalf'; % thresh, splithalf, meanpredict
end
if ~isfield(opts,'elec_max_thresh') || isempty(opts.elec_max_thresh)
    opts.elec_max_thresh = 1; % minimum required maximal response in % signal change for electrode inclusion
end
if ~isfield(opts,'elec_mean_thresh') || isempty(opts.elec_mean_thresh)
    opts.elec_mean_thresh = 0; % minimum required mean response during stim_on period in % signal change
end
if ~isfield(opts,'elec_splithalf_thresh') || isempty(opts.elec_splithalf_thresh)
    opts.elec_splithalf_thresh = -0.2; % minimum required R2 between split halves of data
end
if ~isfield(opts,'elec_meanpredict_thresh') || isempty(opts.elec_meanpredict_thresh)
    opts.elec_meanpredict_thresh = 0; % minimum required R2 for prediction by mean (1 - (SSEresidual/SSEtotal)
end
if ~isfield(opts,'elec_exclude_depth') || isempty(opts.elec_exclude_depth)
    opts.elec_exclude_depth  = false; % boolean
end
if ~isfield(opts,'average_trials') || isempty(opts.average_trials)
    opts.average_trials      = true;  % boolean
end
if ~isfield(opts,'normalize_data') || isempty(opts.normalize_data)
    opts.normalize_data      = false;  % boolean
end
if ~isfield(opts,'average_elecs') || isempty(opts.average_elecs)
    opts.average_elecs       = false; % boolean
end
if ~isfield(opts,'doplots') || isempty(opts.doplots)
    opts.doplots             = true; % boolean
end
if ~isfield(opts,'sort_channels') || isempty(opts.sort_channels)
    opts.sort_channels      = true;  % boolean
end
if ~isfield(opts,'plotsavedir') || isempty(opts.plotsavedir)
    opts.plotsavedir         = 	fullfile(analysisRootPath, 'figures');
end

%% Initialize
savePlots         = opts.doplots;
plotSaveDir       = opts.plotsavedir;
plotSaveDir_epoch = fullfile(plotSaveDir, 'epochselection');
plotSaveDir_elecs = fullfile(plotSaveDir, 'electrodeselection');

% Make figure directories if they don't exist
if savePlots
    if ~exist(plotSaveDir_epoch,'dir'); mkdir(plotSaveDir_epoch); end
    if ~exist(plotSaveDir_elecs,'dir'); mkdir(plotSaveDir_elecs); end
end

% We will use both the voltage and broadband epochs in selection, but will
% output only the broadband epochs for the subsequent model testing step.
allEpochs   = []; % To be filled with broadband epochs across all subjects
allChannels = []; % To be filled with channel info across all subjects

%% Select data

for ii = 1:length(data) % Loop over subjects
    
    subject     = data{ii}.subject;
    epochs_b    = data{ii}.epochs_b;
    epochs_v    = data{ii}.epochs_v;
    t           = data{ii}.t;
    channels    = data{ii}.channels;
    events      = data{ii}.events;

    fprintf('[%s] Selecting data for subject %s \n',mfilename, subject);
          
    % Restrict selection to relevant stimuli only
    stim_idx = contains(events.trial_name, opts.stimnames);
    epochs_b = epochs_b(:, stim_idx, :);
    epochs_v = epochs_v(:, stim_idx, :);
    events = events(stim_idx, :);
    
    % Restrict selection to included channels only
    
    % Exclude depth electrodes
    if opts.elec_exclude_depth 
        chan_idx = ~contains(lower(channels.type), 'seeg');
        channels = channels(chan_idx,:);
        epochs_b = epochs_b(:,:,chan_idx);
        epochs_v = epochs_v(:,:,chan_idx);
    end  
    
    % Include only requested areas
    if ~isempty(opts.areanames)
        if ~iscell(opts.areanames),opts.areanames = {opts.areanames};end
        [chan_groupidx, channels_grouped] = groupElecsByVisualArea(channels);
        chan_idx = [];
        for cc = 1:length(opts.areanames)
            group_idx = strcmp(opts.areanames{cc}, channels_grouped.name);
            chan_idx_matched = find(chan_groupidx{group_idx});
            if ~any(chan_idx_matched)
                fprintf('[%s] Did not find matching electrodes in area %s for subject %s\n', mfilename, opts.areanames{cc}, subject);
            else
                chan_idx = [chan_idx; chan_idx_matched];
            end
        end
        
        channels = channels(chan_idx,:);
        epochs_b = epochs_b(:,:,chan_idx);
        epochs_v = epochs_v(:,:,chan_idx);      

    end
    
    if isempty(channels)
        fprintf('[%s] Did not find any matching electrodes for subject %s. Moving to next subject\n', mfilename, subject);
        continue
    else
        %% STEP 1 Select epochs

        % Exclude bad epochs based on VOLTAGE 
        [epochs_v] = ecog_normalizeEpochs(epochs_v, t, opts.baseline_time, 'subtractwithintrial');

        fprintf('[%s] Removing bad epochs...\n',mfilename);

        %[~, outlier_idx, max_epochs, outlier_thresh] = ecog_selectEpochs(epochs_v, t, opts);
        [~, outlier_idx, max_epochs, outlier_thresh] = ecog_selectEpochsStat(epochs_v, t, opts.stim_on);

        % Plot the included and excluded trials: all channels combined
        if savePlots
            % Make separate plots for different electrode groups
            groups = unique(channels.group);
            nGroups = length(groups);
            for jj = 1:nGroups
                figureName = sprintf('outlierepochs_allchans_sub-%s-%s', subject, groups{jj});
                figure('Name', figureName);hold on
                chan_idx = contains(channels.group, groups{jj});
                outlier_idx_group = outlier_idx(:,chan_idx);
                if any(outlier_idx_group(:)) 
                    subplot(2,2,1); 
                    plot(t, epochs_v(:,outlier_idx_group));
                    title('excluded epochs based on voltage - voltage')
                    subplot(2,2,2); 
                    plot(t, epochs_b(:,outlier_idx_group));
                    title('excluded epochs based on voltage - broadband')
                end
                subplot(2,2,3); 
                plot(t, epochs_v(:,~outlier_idx_group));
                title('included epochs based on voltage - voltage')
                subplot(2,2,4); 
                plot(t, epochs_b(:,~outlier_idx_group));
                title('included epochs based on voltage - broadband')
                % Set axes
                set(findall(gcf,'-property','FontSize'),'FontSize',14)
                set(gcf, 'Position', [150 100 1400 600]);
                % Save
                saveas(gcf, fullfile(plotSaveDir_epoch, figureName), 'png'); close;
            end
        end

        % Plot the outlier trials: individual plots
        if savePlots
            for jj = 1:height(channels)
                % Skip plotting of depth electrodes if those are not included
                if opts.elec_exclude_depth && contains(lower(channels.type(jj)), 'seeg')
                    continue
                end 
                if any(outlier_idx(:,jj))
                    outliers_found = find(outlier_idx(:,jj));
                    nOutliers = length(outliers_found);
                    dim1 = round(sqrt(nOutliers+1)); dim2 = ceil((nOutliers+1)/dim1);
                    dim1 = round((nOutliers+1)/2);
                    dim2 = round((nOutliers+1)/dim1);
                    % voltage
                    figureName = sprintf('outlierepochs_sub-%s_chan-%s-voltage', subject, channels.name{jj});
                    figure('Name', figureName); hold on;
                    subplot(dim2,dim1,1); hold on; title(channels.name{jj}); 
                    histogram(max_epochs(:,jj),100); line([outlier_thresh(jj) outlier_thresh(jj)], get(gca, 'YLim'), 'Color', 'r','LineStyle', ':', 'LineWidth', 2);
                    set(gca, 'fontsize', 14); xlabel('max pows'); ylabel('number of epochs');
                    for kk = 1:nOutliers
                        subplot(dim2,dim1,kk+1); 
                        ecog_plotSingleTimeCourse(t, epochs_v(:,outliers_found(kk),jj), [], [], sprintf('epoch %d %s', outliers_found(kk), events.trial_name{outliers_found(kk)}));    
                    end
                    set(gcf, 'Position', [150 100 300*dim1 300*dim2]);
                    saveas(gcf, fullfile(plotSaveDir_epoch, figureName), 'png'); close;
                    figureName = sprintf('outlierepochs_sub-%s_chan-%s-broadband', subject, channels.name{jj});
                    figure('Name', figureName); hold on;
                    subplot(dim2,dim1,1); hold on; title(channels.name{jj}); 
                    histogram(max_epochs(:,jj),100); line([outlier_thresh(jj) outlier_thresh(jj)], get(gca, 'YLim'), 'Color', 'r','LineStyle', ':', 'LineWidth', 2);
                    set(gca, 'fontsize', 14); xlabel('max pows'); ylabel('number of epochs');
                    for kk = 1:nOutliers
                        subplot(dim2,dim1,kk+1); 
                        ecog_plotSingleTimeCourse(t, epochs_b(:,outliers_found(kk),jj), [], [], sprintf('epoch %d %s', outliers_found(kk), events.trial_name{outliers_found(kk)}));    
                    end
                    set(gcf, 'Position', [150 100 300*dim1 300*dim2]);
                    saveas(gcf, fullfile(plotSaveDir_epoch, figureName), 'png'); close;
                end
            end
        end

        % Mask the broadband epochs to include only the selected epochs.
        epochs = epochs_b;
        epochs(:,outlier_idx) = nan;

        %% STEP 2 Convert to percent signal change 
        fprintf('[%s] Converting epochs to percent signal change...\n',mfilename);

        % Provide run index to perform separately for each run and session
    %     [~,~,task_idx]= unique(events.task_name);
    %     [~,~,ses_idx]= unique(events.session_name);
    %     [~,~,run_idx] = unique(events.run_name);
    %     [~,~,idx] = unique([task_idx ses_idx run_idx], 'rows');
        idx = [];
        [epochs] = ecog_normalizeEpochs(epochs, t, opts.baseline_time, 'percentsignalchange', idx);
        channels.units = repmat({'%change'}, [height(channels),1]);

        %% STEP 3 Select electrodes   
        fprintf('[%s] Selecting electrodes...\n',mfilename);

        [epochs_selected, channels_selected, select_idx, R2, epochs_split] = ...
            ecog_selectElectrodes(epochs, channels, events, t, opts);

        if savePlots     
            if ~isempty(epochs_split)
                for el = 1:height(channels)
                    figureName = sprintf('%s_%s_%s_%s', subject, ...
                    channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});
                    figure;hold on;
                    plot(squeeze(epochs_split(1,el,:)), 'r','LineWidth', 2);
                    plot(squeeze(epochs_split(2,el,:)), 'b','LineWidth', 2);
                    axis tight
                    nSamp = size(epochs,1); nSampTot = nSamp * length(opts.stimnames);
                    set(gca, 'XTick', 1:nSamp:nSampTot, 'XTickLabel', opts.stimnames);
                    xtickangle(45)
                    title(sprintf('%s %s %s R2 = %0.2f', ...
                    channels.name{el}, channels.bensonarea{el}, channels.wangarea{el}, R2(el)));
                    scrSz = get(0, 'Screensize');
                    set(gcf, 'Position', [1 1 scrSz(3)/2 scrSz(4)/2]);
                    set(findall(gcf,'-property','FontSize'),'FontSize',14)
                    saveas(gcf, fullfile(plotSaveDir_elecs, figureName), 'png'); close;
                end
            end

            % Compute mean across all trials
            mean_resp = mean(epochs,2,'omitnan');
            % Compute std deviations for plotting
            llim = (mean_resp - (std(epochs,0,2,'omitnan')));
            ulim = (mean_resp + (std(epochs,0,2,'omitnan')));
            mean_resp_sd = cat(2, llim, ulim);
            % Plot all channels
            nEl = size(mean_resp,3); 
            figureName = sprintf('viselec_%s_all', subject);
            figure('Name', figureName); plotDim1 = round(sqrt(nEl)); plotDim2 = ceil((nEl)/plotDim1);
            for el = 1:nEl
                subplot(plotDim1,plotDim2,el); hold on
                plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        
                ecog_plotSingleTimeCourse(t, mean_resp(:,:,el), squeeze(mean_resp_sd(:,:,el)), [], plotTitle);
                %if el == 1; xlabel('Time (s)'); ylabel('Broadband signal change');end
                set(gcf, 'Position', [150 100 1500 1250]);
            end
            saveas(gcf, fullfile(plotSaveDir_elecs, figureName), 'png'); close;
            % Plot selected channels
            figureName = sprintf('viselec_%s_selected', subject);
            figure('Name', figureName); 
            for el = 1:nEl
                if select_idx(el)
                    subplot(plotDim1,plotDim2,el); hold on
                    plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        
                    ecog_plotSingleTimeCourse(t, mean_resp(:,el), squeeze(mean_resp_sd(:,:,el)), [], plotTitle)    
                    set(gcf, 'Position', [150 100 1500 1250]);
                end
            end
            saveas(gcf, fullfile(plotSaveDir_elecs, figureName), 'png'); close;
        end

        % Select the channels
        epochs = epochs_selected;

        % Update channels table
        channels = channels_selected;

        %% STEP 4 average across trials, concatenate subjects

        % Average trials
        if opts.average_trials
            [epochs] = ecog_averageEpochs(epochs, events, opts.stimnames);   
        end

        % Concatenate the data across subjects
        allEpochs = cat(3, allEpochs, epochs);    

        % Remove a number of columns from channel table for readability, and
        % concatenate across subjects
        channels = removevars(channels, {'low_cutoff', 'high_cutoff', 'reference', 'group', 'bb_method', 'bb_bandwidth'});
        allChannels = [allChannels; channels];
    end   
end

%% Run checks and finalize processing %%

% Now that we have the selected trials and electrodes from all subjects, 
% run a few checks, and make decisions on how to further process the data.
epochs    = allEpochs;
channels  = allChannels;
stimnames = opts.stimnames;

% Check that all channels have same sample frequency
assert(length(unique(channels.sampling_frequency))==1);
srate = channels.sampling_frequency(1);

% Sort electrodes on visual area (rather than subjectID)?
if opts.sort_channels
    % sort on benson area
    [~,I] = sortVisualAreaNames(channels.benson14_varea);
    channels = channels(I,:);
    epochs = epochs(:,:,I);
end

% Scale each electrode to its max?
if opts.normalize_data       
    [epochs] = normalize_data(epochs);
end

% Average elecs within area?
if opts.average_elecs
    fprintf('[%s] Averaging electrodes...\n', mfilename);
    [epochs, channels, epochs_se] = average_elecs(epochs, channels);
else
    epochs_se = [];
end

% Add an index column to channels
index = [1:height(channels)]';
channels = addvars(channels, index, 'Before', 'name'); 

fprintf('[%s] Done! \n',mfilename);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
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
function [data, channels, se] = average_elecs(data, channels)
    
    area_names = {'V1', 'V2', 'V3', 'V3ab', 'LOTO', 'IPS'};
    [chan_idx, channels, group_prob] = groupElecsByVisualArea(channels, 'probabilisticresample', area_names);
    fun = @mean;
    [data, se] = averageWithinArea(data, group_prob, fun);   
end
