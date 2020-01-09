% s_determineOnsetShiftUMCUvsNYU
%
% Purpose of this script:
%
% Visual comparison of the broadband responses from V1 electrodes in the
% UMCU patient (sub-chaam) and NYU patients suggested a systematic offset
% in the onset of the visual response between the two sites.
%
% We concluded this is due to a delay in the stimulus presentation relative
% to the trigger at UMCU due to differences in the experimental setup. We
% decide to shift the stimulus onsets for the UMCU data by a fixed amount
% to bring the two datasets in alignment.
%
% The shift is determined through cross-correlation of the evoked responses
% between sub-chaam (UMCU) versus sub-som708 and sub-692 (NYU).
%
% We perform the following steps:
% 1. Loading in the common-averaged referenced data (voltage timecourses)
%   from the BIDS derivatives folder for each subject.
% 2. Use the benson atlas to selecting electrodes that are within V1/V2 and
%   that have a predicted eccentricity of less than 10 degrees. 
% 3. Computing the ERPs (average response across all trials), performing
%   baseline correction using the prestimulus period.
% 4. Compute cross correlation between all electrodes of sub-chaam and the
%   two NYU subjects, separately for V1 and V2. 
% 5. Compute the average shift that is needed to maximize the correlation
%   between UMCU and NYU electrodes, separately for V1 and V2.
% 6. Take the average across V1 and V2 as the estimated onset shift.
%
% Iris Groen, Dec 2019, BAIR.


%% Set up variables

subjectList             = {'chaam', 'som708', 'som692'}; 
inputDir                = fullfile(bidsRootPath, 'derivatives', 'ECoGCAR'); % where the evoked data are
tasks                   = {'spatialpattern', 'temporalpattern', 'soc'};
description             = 'reref';
epochTime               = [-0.1 0.5];
baselineTime            = [-0.1 0];
eccentricityThreshold   = 10;

%% Load the data from patients with V1 coverage
data = cell(length(subjectList),1);

for ii = 1 : length(subjectList)
    
    subject = subjectList{ii};
	fprintf('[%s] Processing subject %s...\n',mfilename, subject);

    % Get visual area matches for this subject
    fprintf('[%s] Computing matches with visual atlases...\n',mfilename);
    specs = [];
    specs.pID           = subject; 
    specs.plotmesh      = 'none';
    BIDSformatted       = 1;
    visualelectrodes    = electrode_to_nearest_node(specs, BIDSformatted);

    % Get the data
    fprintf('[%s] Loading data...\n',mfilename);
    [subdata, channels, events] = bidsEcogGetPreprocData(inputDir, subject, [], tasks, [], description);
    
    % Select a subset of channels
    fprintf('[%s] Selecting electrodes...\n',mfilename);

    % Add visual area names (W and B) ecc, angle, sigma to channels table
    [channels] = bair_addVisualAtlasNamesToChannelTable(channels,visualelectrodes);
    
    chan_idx = find(contains(channels.bensonarea, {'V1', 'V2'}) & contains(channels.status, 'good') & channels.bensoneccen < eccentricityThreshold);                
    
    subdata = subdata(chan_idx,:);
	channels = channels(chan_idx,:);
        
    % Resample the UMCU data 
    if contains(subject, {'chaam'})
        fprintf('[%s] This is a umcu patient: resampling \n',mfilename);

        % Resample data; assuming desired sample rate of 512
        subdata = downsample(subdata', channels.sampling_frequency(1)/512)';
        events.event_sample = round(events.event_sample/(channels.sampling_frequency(1)/512));
        channels.sampling_frequency(:) = 512;
        
        % Remove electrodes identified as epileptic % TEMPORARY UNTIL GIO
        % FIXES BIDS FORMATTING FOR THIS PATIENT AFTER WHICH THESE
        % ELECTRODES SHOULD BE INDICATED AS "BAD" IN THE ELECTRODES FILES
        chan_idx = find(~contains(channels.name, {'Oc12', 'Oc13', 'Oc21', 'Oc22'}));
        subdata = subdata(chan_idx,:);
        channels = channels(chan_idx,:);
    end

    % Epoch the data
    fprintf('[%s] Epoching data...\n',mfilename);
    [epochs, t] = ecog_makeEpochs(subdata, events.event_sample, epochTime, channels.sampling_frequency(1));  

    fprintf('[%s] Found %d epochs across %d runs and %d sessions \n', ...
        mfilename, size(epochs,2), length(unique(events.run_name)), length(unique(events.session_name)));

    % Perform baseline correction
    fprintf('[%s] Baseline correcting epochs...\n',mfilename);
    [epochs] = ecog_normalizeEpochs(epochs, t, baselineTime, 'subtractwithintrial');
    channels.units = repmat({'%change'}, [height(channels),1]);
    
    % Add a subject index column to events and channels
    events.subject_name = repmat({subject}, [height(events),1]);
    channels.subject_name = repmat({subject}, [height(channels),1]);
        
    % Collect in struct
    data{ii}.subject  = subject;
	data{ii}.epochs   = epochs;
    data{ii}.t        = t;
    data{ii}.events   = events;
    data{ii}.channels = channels;

end

%% Compute average timecourses to compare
avData = [];
avChanInfo = [];
chanTypes = {'V1', 'V2'};

for ii = 1 : length(subjectList)
    for jj = 1 : length(chanTypes)
        chanInx = find(contains(data{ii}.channels.bensonarea, chanTypes{jj}));
        if ~isempty(chanInx)
            dataToAdd = squeeze(mean(data{ii}.epochs(:,:,chanInx),2));
            dataToAdd = dataToAdd ./ max(abs(dataToAdd),[],1);
            avData = [avData dataToAdd];
            avChanInfo = [avChanInfo; [data{ii}.channels.name(chanInx,:) data{ii}.channels.bensonarea(chanInx,:) data{ii}.channels.subject_name(chanInx,:)]];
        end
    end
end

avChanInfo = cell2table(avChanInfo, 'VariableNames', {'channelname', 'bensonarea', 'subjectname'});


% Flip the polarity for patient chaam so that the first peak is negative
% for all channels (just for visualization, cross-correlations are computed
% on absolute values; see below)
for ii = 1 : size(avData,2)
    if contains(avChanInfo.subjectname{ii}, 'chaam') & ...
            contains(avChanInfo.channelname{ii}, {'Oc12', 'Oc28', 'Oc17', 'Oc27', 'Oc26'})
        avData(:, ii) = -1 * avData(:,ii);
    end
end

%% Plot the data
t = data{1}.t;

figure('Name', 'Selected data'); 
for ii = 1 : length(subjectList)
    subplot(1,length(subjectList),ii); hold on
    l = []; 
    for jj = 1 :length(chanTypes)
        chanInx = find(contains(avChanInfo.subjectname, subjectList{ii}) & contains(avChanInfo.bensonarea, chanTypes{jj}));
        if ~isempty(chanInx)
            if jj == 1, colors = autumn(length(chanInx)); else, colors = winter(length(chanInx)); end
            dataToPlot = avData(:,chanInx);
            for kk = 1:size(dataToPlot,2)
                plot(t,dataToPlot(:,kk), 'LineWidth', 2, 'Color', colors(kk,:));
                l = [l; {sprintf('%s - %s', avChanInfo.channelname{chanInx(kk)}, avChanInfo.bensonarea{chanInx(kk)})}];
            end 
        end
    end
    yaxlims = get(gca, 'YLim');
    line([0 0], [yaxlims(1) yaxlims(2)], 'Color', 'k', 'LineStyle', ':')
    line([t(1) t(end)], [0 0],'Color', 'k', 'LineStyle', ':');
    legend(l)
	title(subjectList{ii});
    xlabel('time');ylabel('ERP amplitude');
    set(gca,'FontSize', 14);
end
set(gcf, 'Position', [800 750 1750 600]);

%% Compute cross correlations between UMCU and NYU for V1 and V2

lagDiff = []; l = [];
figure('Name', 'Cross-correlations'); 
for ii = 1 : length(chanTypes)
    subplot(2,1,ii); title(chanTypes{ii});  hold on;
    set(gcf, 'Position', [0 200 800 1200]);
    set(gca,'FontSize', 14);
    c = 1;
    chanInx1 = find(contains(avChanInfo.bensonarea, chanTypes{ii}) & contains(avChanInfo.subjectname, 'chaam'));
    chanInx2 = find(contains(avChanInfo.bensonarea, chanTypes{ii}) & ~contains(avChanInfo.subjectname, 'chaam'));
    for jj = 1 : length(chanInx1)
        for kk = 1 : length(chanInx2)
            % Compute cross correlation between the absolute time courses
            s1 = abs(avData(:,chanInx1(jj)));
            s2 = abs(avData(:,chanInx2(kk)));
            [acor,lag] = xcorr(s1,s2);
            % Plot
            plot(lag, acor, 'LineWidth', 2);
            [~,I] = max(abs(acor));
            lagDiff{ii}(c) = lag(I);
            comparisonName = sprintf('%s %s - %s %s', avChanInfo.subjectname{chanInx1(jj)}, avChanInfo.channelname{chanInx1(jj)}, ...
                                            avChanInfo.subjectname{chanInx2(kk)}, avChanInfo.channelname{chanInx2(kk)});
            disp(comparisonName);                            
            disp((lag(I)/512)*1000);
            l{ii}{c} = comparisonName;
            c = c + 1;
            legend(l{ii});
        end        
    end
    xlabel('lag (samples)'); ylabel('cross-correlation');
end

%% Compute shift

shiftInSamples = [];
shiftInSeconds = [];
for ii = 1 : length(chanTypes)
    shiftInSamples(ii) = median(lagDiff{ii})+1; % Add one to account to the lag difference to account for 1-based indexing in Matlab
    shiftInSeconds(ii) = shiftInSamples(ii)/512;
	fprintf('Estimated shift for %s = %d samples, %0.3f seconds \n', chanTypes{ii}, round(shiftInSamples(ii)), shiftInSeconds(ii));
end

% Take mean across V1 and V2
shiftInSamples = mean(shiftInSamples);
shiftInSeconds = shiftInSamples/512;
fprintf('Mean estimated shift = %d samples, %0.3f seconds \n', round(shiftInSamples), shiftInSeconds);

%% Plot shift

figure('Name', 'Selected data - shifted'); 
set(gcf, 'Position', [800 750 1500 800]);
colors = jet(length(subjectList));

% Original data
for ii = 1 : length(chanTypes)
    subplot(2,2,ii); title(chanTypes{ii}); hold on;
    set(gca, 'FontSize', 14);
    for jj = 1 : length(subjectList)
        chanInx = find(contains(avChanInfo.bensonarea, chanTypes{ii}) & contains(avChanInfo.subjectname, subjectList{jj}));
        if ~isempty(chanInx)
            p = plot(t,avData(:,chanInx), 'LineWidth', 2, 'Color', colors(jj,:), 'DisplayName', subjectList{jj});
        end
        legend;
    end
    yaxlims = get(gca, 'YLim');
    l = line([0 0], [yaxlims(1) yaxlims(2)], 'Color', 'k', 'LineStyle', ':');
    l.Annotation.LegendInformation.IconDisplayStyle = 'off';
    l = line([t(1) t(end)], [0 0],'Color', 'k', 'LineStyle', ':');
    l.Annotation.LegendInformation.IconDisplayStyle = 'off';
    xlabel('Time (ms)'); ylabel('ERP amplitude (normalized)');

end

% Shifted data
for ii = 1 : length(chanTypes)
    subplot(2,2,ii+2); title([chanTypes{ii} ' shifted by ' sprintf('%2.0f ms', shiftInSeconds*1000)]); hold on;
    set(gca, 'FontSize', 14);
    
    for jj = 1 : length(subjectList)
        chanInx = find(contains(avChanInfo.bensonarea, chanTypes{ii}) & contains(avChanInfo.subjectname, subjectList{jj}));
        if ~isempty(chanInx)
            if contains(subjectList{jj}, 'chaam')
                timeToPlot = t-shiftInSeconds;
            else
                timeToPlot = t;
            end
            p = plot(timeToPlot,avData(:,chanInx), 'LineWidth', 2, 'Color', colors(jj,:), 'DisplayName', subjectList{jj});
        end
        legend;
    end
    yaxlims = get(gca, 'YLim');
    l = line([0 0], [yaxlims(1) yaxlims(2)], 'Color', 'k', 'LineStyle', ':');
    l.Annotation.LegendInformation.IconDisplayStyle = 'off';
    l = line([t(1) t(end)], [0 0],'Color', 'k', 'LineStyle', ':');
    l.Annotation.LegendInformation.IconDisplayStyle = 'off';
    set(gca, 'XLim', epochTime);
    xlabel('Time (ms)'); ylabel('ERP amplitude (normalized)');
end


