
% Check whether we have the ECoG_utils repository on the path
if ~exist('ecog_matchChannels', 'file')
    tbUse ECoG_utils;
end

addpath(genpath('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal'));
addpath(genpath('/Volumes/server/Projects/Temporal_integration/DN_2018_data/code'))

stimNm = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5',...
'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};

saveFigure = 1;
%figureFormat = 'PNG';  % PNG or EPS

saveLoc = fullfile(dn_ctrst_RootPath, 'dataFigures');
    
%% load the data

fileName = 'preproc_epoched_visualelecs.mat';
load(fullfile('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess', fileName));

%% compute baseline a la Zhou et al 2019
t = a{1}.trials.time;

% CHANGE BROADBAND DATA BASELINE
bb = {};
for k = 1 : length(a)
    bb{k} = a{k}.trials.broadband;
end

% DEFINE BASELINE RANGE ---------------------------------------------------
% use the firsrt 200ms before stimulus onset as the baseline
base_range = (t >= -0.2 & t < 0);

for k = 1 : length(a)
    m_base{k} = squeeze(median(mean(bb{k}(:, base_range, :), 2), 3));
    bb{k}     = bb{k}./m_base{k}-1;
end

%% select elecs based on criteria Zhou et al 2019

selTresh = 1;
stim_on = t > 0 & t < 0.2;
selectIdx = cell(size(a));

for k = 1 : length(a)
    
    % use contrast temporal stimuli
    stimsForSelection = contains(a{k}.trials.events.trial_name, stimNm);
    
    % use HRF task, if not available, use other trials from spatialpattern

%     stimsForSelection = contains(a{k}.trials.events.trial_name, {'pattern', 'HRFPATTERN', 'HRF', 'hrf'});
%     if max(stimsForSelection) == 0
%         % this patient has no HRF; use other trials from spatialpattern
%         stimsForSelection = contains(a{k}.trials.events.trial_name, {'SPARSITY', 'PLAID', 'CIRCULAR'});
%         stim_on = t > 0 & t < 0.2;
%     end

    mbb = mean(bb{k}(:,:,stimsForSelection),3);
    llim = (mbb - (std(bb{k}(:,:,stimsForSelection),0,3)));
    ulim = (mbb + (std(bb{k}(:,:,stimsForSelection),0,3)));
    mbb_ci = cat(3, llim, ulim);
    nEl = size(mbb,1); 
    
    % plot
    figureName = sprintf('viselec_%s_%s_all', a{k}.sub, a{k}.ses);
    figure('Name', figureName); nSubPlot = ceil(sqrt(nEl)); 
    for el = 1:nEl
        subplot(nSubPlot,nSubPlot,el); hold on
        plotTitle = sprintf('%s W:%s B:%s ', a{k}.trials.channels.name{el}, a{k}.trials.channels.wangatlas{el}, a{k}.trials.channels.bensonatlas{el});        
        ecog_plotSingleTimeCourse(t, mbb(el,:), squeeze(mbb_ci(el,:,:)), [], plotTitle)       
        set(gcf, 'Position', [150 100 1500 1250]);
    end
    if saveFigure
        %printnice([], 0, saveLoc, figureName)
        saveas(gcf, fullfile(saveLoc, 'electrode_selection', figureName), 'png'); close;
    end
    
    % FIRST SELECTION CRITERION:
    select_idx1 = max(mbb(:,stim_on),[],2) > selTresh;
    
    % SECOND SELECTION CRITERION:
    select_idx2 = zeros(size(select_idx1));
    for el = 1:nEl
        if mean(mbb(el, stim_on)) > 0
            select_idx2(el) = 1;
        else
            select_idx2(el) = 0;
        end
    end
    
    % EXCLUDE CHANNELS LABELED AS BAD
    if isfield(summary(a{k}.trials.channels), 'status')
        select_idx3 = contains(a{k}.trials.channels.status, 'good');
    else
        select_idx3 = ones(size(select_idx1));
    end
    
    % EXCLUDE DEPTH ELECTRODES
    select_idx4 = contains(a{k}.trials.channels.type, {'ecog', 'ECOG'});
    
    % combine criteria
    select_idx = select_idx1 + select_idx2 + select_idx3 + select_idx4;
    select_idx = (select_idx == 4);
    
    % plot selection
    figureName = sprintf('viselec_%s_%s_selected', a{k}.sub, a{k}.ses);
    figure('Name', figureName); nSubPlot = ceil(sqrt(nEl)); 
    for el = 1:nEl
        if select_idx(el)
            subplot(nSubPlot,nSubPlot,el); hold on
            plotTitle = sprintf('%s W:%s B:%s ', a{k}.trials.channels.name{el}, a{k}.trials.channels.wangatlas{el}, a{k}.trials.channels.bensonatlas{el});
            ecog_plotSingleTimeCourse(t, mbb(el,:), squeeze(mbb_ci(el,:,:)), [], plotTitle)    
            set(gcf, 'Position', [150 100 1500 1250]);
        end
    end
    if saveFigure
        %printnice([], 0, saveLoc, figureName)
        saveas(gcf, fullfile(saveLoc, 'electrode_selection', figureName), 'png'); close;
    end    
    selectIdx{k} = select_idx;
end

%% average trials across condition, and group by visual area
    
data = [];
se = [];

% concatenate all selected electrodes
for k = 1:length(a)
    temp_data = [];
    temp_se = [];
    for k1 = 1:length(stimNm)
        trialIdx = contains(a{k}.trials.events.trial_name, stimNm{k1});
        %temp(:,k1,:) = squeeze(median(bb{k}(selectIdx{k},:,trialIdx), 3));
        temp_data(:,k1,:) = squeeze(mean(bb{k}(selectIdx{k},:,trialIdx), 3));
        temp_se(:,k1,:) = squeeze(std(bb{k}(selectIdx{k},:,trialIdx),0, 3))/sqrt(length(find(trialIdx)));
    end
    data = cat(1, data, temp_data);
	se = cat(1, se, temp_se);
end

% collect info about elecs
patientname = {};
sessionname = {};
elecname = {};
wang = {};
benson = {};

for k = 1:length(a)
    idx = length(patientname)+1:length(patientname)+length(find(selectIdx{k}));
    if ~isempty(idx)
        patientname(idx,1) = repmat({a{k}.sub}, [length(idx) 1]);
        sessionname(idx,1) = repmat({a{k}.ses}, [length(idx) 1]);
        elecname(idx,1) = a{k}.trials.channels.name(selectIdx{k});
        wang(idx,1) = a{k}.trials.channels.wangatlas(selectIdx{k});
        benson(idx,1) = a{k}.trials.channels.bensonatlas(selectIdx{k});
    end
end

elec_info = table(patientname, sessionname, elecname, wang, benson);

%% VISUALIZE DATA

for ee = 1:size(data,1)
    figureName = sprintf('contrasttemporal_%s_%s_%s_Wang%s_Benson%s', elec_info.patientname{ee}, elec_info.sessionname{ee}, elec_info.elecname{ee}, elec_info.wang{ee}, elec_info.benson{ee});
    figure('Name', figureName);
    %y_limits = [-1 10];
    y_limits = [-1 max(data(ee,:))];
	m = squeeze(data(ee,:,:));
    s = squeeze(se(ee,:,:));
    for ii = 1:5 % CRF
       subplot(3,6,ii)
       ecog_plotSingleTimeCourse(t, m(ii,:), [m(ii,:)-s(ii,:);  m(ii,:)+s(ii,:)],[],stimNm{ii},[],y_limits);   
    end
    for ii = 6:11
       subplot(3,6,ii+1)
       ecog_plotSingleTimeCourse(t, m(ii,:), [m(ii,:)-s(ii,:);  m(ii,:)+s(ii,:)],[],stimNm{ii},[],y_limits);   
    end
    for ii = 12:17
       subplot(3,6,ii+1)
       ecog_plotSingleTimeCourse(t, m(ii,:), [m(ii,:)-s(ii,:);  m(ii,:)+s(ii,:)],[],stimNm{ii},[],y_limits);   
    end
    set(gcf, 'Position', [150 100 1500 1250]);
    if saveFigure
        %printnice([], 0, saveLoc, figureName)
        saveas(gcf, fullfile(saveLoc, 'selectedelectrodes_responsebycondition', figureName), 'png'); close;
    end
end

%% average by areas
AreaNames = {'V1', 'V2', 'V3', 'Vhigher'};

INX = [];
INX{1} = contains(elec_info.wang, 'V1') | contains(elec_info.benson, 'V1');
INX{2} = contains(elec_info.wang, 'V2') | contains(elec_info.benson, 'V2');
INX{3} = contains(elec_info.wang, 'V3') & ~contains(elec_info.wang, {'V3a', 'V3b'}) | contains(elec_info.benson, 'V3') & ~contains(elec_info.benson, {'V3a', 'V3b'});
INX{4} = contains(elec_info.wang, {'V3a', 'V3b', 'IPS', 'hV4', 'LO', 'TO'}) | contains(elec_info.benson, {'V3a', 'V3b', 'hV4', 'LO', 'TO'}) ;

for vv = 1:1%length(AreaNames)
    figureName = sprintf('averagebycondition_%s',AreaNames{vv});
    figure('Name', figureName);
    ee = INX{vv};
    mdata = squeeze(mean(data(ee,:,:),1));
    sdata = squeeze(std(data(ee,:,:),0,1))/sqrt(length(find(ee)));
    %y_limits = [-1 max(mdata(:))];
    y_limits = [-1 25];
    for ii = 1:5 % CRF
       subplot(3,6,ii)
       ecog_plotSingleTimeCourse(t, mdata(ii,:), [mdata(ii,:) - sdata(ii,:); mdata(ii,:) + sdata(ii,:)],[],stimNm{ii},[],y_limits);   
    end
    for ii = 6:11 % ONEPULSE
       subplot(3,6,ii+1)
       ecog_plotSingleTimeCourse(t, mdata(ii,:), [mdata(ii,:) - sdata(ii,:); mdata(ii,:) + sdata(ii,:)],[],stimNm{ii},[],y_limits);   
    end
    for ii = 12:17 % TWOPULSE
       subplot(3,6,ii+1)
       ecog_plotSingleTimeCourse(t, mdata(ii,:), [mdata(ii,:) - sdata(ii,:); mdata(ii,:) + sdata(ii,:)],[],stimNm{ii},[],y_limits);   
    end
    set(gcf, 'Position', [150 100 1500 1250]);
    if saveFigure
        %printnice([], 0, saveLoc, figureName)
        %saveas(gcf, fullfile(saveLoc, 'selectedelectrodes_responsebycondition', figureName), 'png'); close;
    end
end

%% save
saveName = 'data_contrasttemporal_thresh1.mat';
save(fullfile('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess', saveName), 'data', 'elec_info','t', '-v7.3');

%% save figures for OHBM poster

load(fullfile('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess', saveName));

% add stim boxcar behind plot
contrasts = [0.0625 0.125 0.25 0.5 1]; 
durs      = [0.016667, 0.033333, 0.066667, 0.13333, 0.26667, 0.53333];

nStim = size(data, 2);
stim  = zeros(nStim, length(t));

% CONTRAST STIMULI --------------------------------------------------------
stim(1 : 5, t > 0 & t<=0.5) = 1;
for k = 1 : 5, stim(k, :) = stim(k, :) .* contrasts(k); end

% INCREASING DURATIONS ----------------------------------------------------
for k = 1 : 6, stim(k + 5, (t>0) & (t <= durs(k))) = 1; end

% INCREASING ISI ----------------------------------------------------------
stim(12 : nStim, t > 0 & t <= durs(4)) = 1;
for k = 1 : 6
    t_start = durs(4) + durs(k);
    t_end   = durs(4) * 2 + durs(k);
    stim(11 + k, t > t_start & t <= t_end) = 1;
end

% % VISUALIZE THE STIMULI ---------------------------------------------------
% figure (2)
% for k = 1 : 17
%    subplot(17, 1, k)
%    plot(t, stim(k, :), 'k-')
% end
stim = stim*25;
for vv = 1:1
    figureName = sprintf('averagebycondition_%s',AreaNames{vv});
    figure('Name', figureName);
    ee = INX{vv};
    mdata = squeeze(mean(data(ee,:,:),1));
    sdata = squeeze(std(data(ee,:,:),0,1))/sqrt(length(find(ee)));
    %y_limits = [-1 max(mdata(:))];
    y_limits = [-1 25];
    colors = cool(6);
    for ii = 6:11 % ONEPULSE
       subplot(3,6,ii-5)
      
       %figureName = sprintf('trialaverage_%s_%s',AreaNames{vv}, stimNm{ii});
       %figure('Name', figureName);
       ecog_plotSingleTimeCourse(t, mdata(ii,:), [mdata(ii,:) - sdata(ii,:); mdata(ii,:) + sdata(ii,:)],colors(ii-5,:),[],[],y_limits);
       hold on;
       s = plot(t,stim(ii,:), 'k-'); s.Annotation.LegendInformation.IconDisplayStyle = 'off';
        if ii == 6
            xlabel('Time (s)'); 
            ylabel('Increase in broadband (xfold)');
        else
            set(gca, 'XTick', [], 'YTick', []);
       end
       %printnice([], 0, fullfile(saveLoc, 'selectedelectrodes_responsebycondition'), figureName);close
    end
    for ii = 12:17 % TWOPULSE
       subplot(3,6,ii-5)
       %figureName = sprintf('trialaverage_%s_%s',AreaNames{vv}, stimNm{ii});
       %figure('Name', figureName);
       ecog_plotSingleTimeCourse(t, mdata(ii,:), [mdata(ii,:) - sdata(ii,:); mdata(ii,:) + sdata(ii,:)],colors(ii-11,:),[],[],y_limits);  
       hold on;
       s = plot(t,stim(ii,:), 'k-'); s.Annotation.LegendInformation.IconDisplayStyle = 'off';
       set(gca, 'XTick', [], 'YTick', []);
       %printnice([], 0, fullfile(saveLoc, 'selectedelectrodes_responsebycondition'), figureName);close
    end
    colors = cool(5);

    for ii = 1:5 % CRF
       subplot(3,6,ii+12)
       %figureName = sprintf('trialaverage_%s_%s',AreaNames{vv}, stimNm{ii});
       %figure('Name', figureName);
       ecog_plotSingleTimeCourse(t, mdata(ii,:), [mdata(ii,:) - sdata(ii,:); mdata(ii,:) + sdata(ii,:)],colors(ii,:),[],[],y_limits); 
       hold on;
       s = plot(t,stim(ii,:), 'k-'); s.Annotation.LegendInformation.IconDisplayStyle = 'off';
        set(gca, 'XTick', [], 'YTick', []);
       %printnice([], 0, fullfile(saveLoc, 'selectedelectrodes_responsebycondition'), figureName);close
    end
end
set(gcf, 'Position', [0 0 2500 1250]);

printnice([], 0, saveLoc, figureName)


