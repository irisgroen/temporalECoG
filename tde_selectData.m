function [epochs, channels, stimNames, t, srate] = tde_selectData(data, stimNames, opts)

% Description
%
% function [epochs, channels, stimNames, t] = tde_selectData(data, [stimNames], [opts])
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
if ~exist('data', 'var') || isempty(data)
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

% Initialize
allEpochs   = [];
allChannels = [];

%% Loop over subjects

for ii = 1:length(data) % Loop over subjects
    
    subject     = data{ii}.subject;
    epochs      = data{ii}.epochs;
    channels    = data{ii}.channels;
    events      = data{ii}.events;
    t           = data{ii}.t;

    fprintf('[%s] Selecting data for subject %s \n',mfilename, subject);
          
    % Restrict selection to relevant stimuli only
    stimsForSelection = contains(events.trial_name, stimNames);
    epochs = epochs(:, stimsForSelection, :);
    events = events(stimsForSelection, :);
    
%% STEP 1 Select epochs
    
    %fprintf('[%s] Removing bad epochs...\n',mfilename);
    [epochs_selected, outliers, max_epochs] = ecog_selectEpochs(epochs, t, opts.stim_on, opts.epoch_outlier_thresh);
    
    if savePlots
        for jj = 1:height(channels)
            if any(outliers(:,jj))
                figureName = sprintf('outlierepochs_sub-%s_chan-%s', subject, channels.name{jj});
                figure('Name', figureName); hold on;
                outliers_found = find(outliers(:,jj));
                nOutliers = length(outliers_found);
                dim1 = round((nOutliers+1)/2);
                dim2 = round((nOutliers+1)/dim1);
                subplot(dim2,dim1,1); hold on; title(channels.name{jj}); 
                histogram(max_epochs(:,jj),100); line([opts.epoch_outlier_thresh opts.epoch_outlier_thresh], get(gca, 'YLim'), 'Color', 'r','LineStyle', ':', 'LineWidth', 2);
                set(gca, 'fontsize', 14); xlabel('max broadband'); ylabel('number of epochs');
                for kk = 1:nOutliers
                    subplot(dim2,dim1,kk+1); 
                    ecog_plotSingleTimeCourse(t, epochs(:,outliers_found(kk),jj), [], [], sprintf('epoch %d %s', outliers_found(kk), events.trial_name{outliers_found(kk)}));    
                end
                set(gcf, 'Position', [150 100 300*dim1 300*dim2]);
                saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
            end
        end
    end
    
    epochs = epochs_selected;
    
%% STEP 2 Convert to percent signal change 
    %fprintf('[%s] Converting epochs to percent signal change...\n',mfilename);

    % Provide run index to perform separately for each run and session
%     [~,~,task_idx]= unique(events.task_name);
%     [~,~,ses_idx]= unique(events.session_name);
%     [~,~,run_idx] = unique(events.run_name);
%     [~,~,idx] = unique([task_idx ses_idx run_idx], 'rows');
    idx = [];
    [epochs] = ecog_normalizeEpochs(epochs, t, opts.baseline_time, 'percentsignalchange', idx);
    channels.units = repmat({'%change'}, [height(channels),1]);
    
%% STEP 3 Select electrodes   
    %fprintf('[%s] Selecting electrodes...\n',mfilename);
    
    % Compute mean across all trials
    mean_resp = mean(epochs,2, 'omitnan');
    
    % Initialize selection to include all channels
    select_idx = ones(height(channels),3);
    
    % Exclude channels based on max and mean:
    stim_on_idx = t > opts.stim_on(1) & t <= opts.stim_on(2);  
    select_idx(:,1) = max(mean_resp(stim_on_idx,:),[],1) > opts.elec_max_thresh;
    select_idx(:,2) = mean(mean_resp(stim_on_idx,:),1) > opts.elec_mean_thresh;

    % Exclude depth electrodes
    if opts.elec_exclude_depth, select_idx(:,3) = contains(lower(channels_in.type), 'ecog'); end  

    % Combine criteria
    select_idx = sum(select_idx,2) == 3;

    if savePlots
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
            plotTitle = sprintf('%s W:%s B:%s ', channels.name{el}, channels.wangarea{el}, channels.bensonarea{el});        
            ecog_plotSingleTimeCourse(t, mean_resp(:,:,el), squeeze(mean_resp_sd(:,:,el)), [], plotTitle);
            %if el == 1; xlabel('Time (s)'); ylabel('Broadband signal change');end
            set(gcf, 'Position', [150 100 1500 1250]);
        end
        saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
        % Plot selected channels
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
    
    % Average trials
    if opts.average_trials
        [epochs] = average_trials(epochs, events, stimNames);
    end
    
    % Concatenate the data across subjects
    allEpochs = cat(3, allEpochs, epochs);    
    
    % Remove a number of columns from channel table for readability, and
    % concatenate across subjects
    channels = removevars(channels, {'low_cutoff', 'high_cutoff', 'reference', 'group', 'bb_method', 'bb_bandwidth'});
    allChannels = [allChannels; channels];
    
end

%% Run checks and finalize processing %%

% Now that we have the selected trials and electrodes from all subjects, 
% run a few checks, and make decisions on how to further process the data.
epochs   = allEpochs;
channels = allChannels;

% Check that all channels have same sample frequency
assert(length(unique(channels.sampling_frequency))==1);
srate = channels.sampling_frequency(1);

% Sort electrodes on visual area (rather than subjectID)?
if opts.sort_channels
    [channels] = sort_channels(channels);
end

% Scale each electrode to its max?
if opts.normalize_data       
    [epochs] = normalize_data(epochs);
end

% Average elecs within area?
if opts.average_elecs
    [epochs, channels] = average_elecs(epochs, channels);
end

% Add an index column to channels
index = [1:height(channels)]';
channels = addvars(channels, index, 'Before', 'name'); 

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
        epochs_averaged(:,ii,:) = mean(epochs(:,trial_idx,:),2, 'omitnan');
    end
    % Also compute se across trials?
end

% Channel sorting
function [channels] = sort_channels(channels)   
    % sort on benson area
    [~,I] = sortVisualAreaNames(channels.bensonarea);
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
    subjects = {}; nelecs = {};
    
    INX = [];
    INX{1} = contains(channels.wangarea, 'V1') | contains(channels.bensonarea, 'V1');
    INX{2} = contains(channels.wangarea, 'V2') | contains(channels.bensonarea, 'V2');
    INX{3} = contains(channels.wangarea, 'V3') & ~contains(channels.bensonarea, {'V3a', 'V3b'}) | contains(channels.bensonarea, 'V3') & ~contains(channels.bensonarea, {'V3a', 'V3b'});
    INX{4} = contains(channels.wangarea, {'V3a', 'V3b'}) | contains(channels.bensonarea, {'V3a', 'V3b'}) ;
    INX{5} = contains(channels.wangarea, {'hV4'}) | contains(channels.bensonarea, {'hV4'});
    INX{6} = contains(channels.wangarea, {'LO1', 'LO2'}) | contains(channels.bensonarea, {'LO1', 'LO2'});
    INX{7} = contains(channels.wangarea, {'TO1', 'TO2'}) | contains(channels.bensonarea, {'TO1', 'TO2'});
    INX{8} = contains(channels.wangarea, {'IPS0', 'IPS1', 'IPS2', 'IPS3', 'IPS4', 'IPS5'});
    INX{9} = contains(channels.wangarea, {'VO1','VO2', 'PHC1', 'PHC2'}) | contains(channels.bensonarea, {'VO1', 'VO2'});

    for ii = 1:length(INX)
        mdata = mean(data(:,:,INX{ii}),3);
        
        if opts.normalize_data     
            % another normalization step from Jings code, necessary?
            maxmdata = max(mdata(:));
            mdata    = mdata ./maxmdata;
        end
        
        avdata(:,:,ii) = mdata;
        temp = unique(channels.subject_name(INX{ii}));
        subjects{ii} = [temp{:}];
        nelecs{ii} = length(find(INX{ii}));
    end
    data = avdata;
    
    % Create a new channels table:
    name               = {'V1', 'V2', 'V3', 'V3a V3b', 'hV4', 'LO','TO' 'IPS', 'VO PHC'}';
    type               = repmat({'n/a'}, [length(name) 1]);
    units              = repmat(channels.units(1), [length(name) 1]);
    sampling_frequency = repmat(channels.sampling_frequency(1), [length(name) 1]);
    subject_name       = subjects';    
    number_of_elecs    = nelecs';
    channels = table(name, type, units, sampling_frequency, subject_name, number_of_elecs);
    
end
end