function [out] = tde_selectData(data, stimNames, baselineTime, saveFigure)
% Description      


% Loop over subjects
for ii = 1:length(data)
    epochs   = data.epochs;
    channels = data.channels;
    events   = data.events;

    %% STEP 1 CHECK FOR OUTLIERS EPOCHS/CHANNELS
        % fprintf('[%s] Checking for bad epochs ...\n',mfilename);

        % compute sum over time for each epoch?
        % remove epochs that have 20x the average
        % if more than X epochs for a channel, remove entire channel
        % output a description of how many trials were removed (write to
        % file?)
        
  %% STEP 2 convert to percent signal change 
        fprintf('[%s] Converting epochs to percent signal change...\n',mfilename);

        % Provide run index to perform separately for each run and session
        [~,~,ses_idx]= unique(events.session_name);
        [~,~,run_idx] = unique(events.run_name);
        idx = (ses_idx*100)+run_idx;
        
        [epochs] = ecog_normalizeEpochs(epochs, t, baselineTime, 'percentsignalchange', idx);
         
 %% STEP 3 select electrodes
        
%         selTresh = 1;
%         stim_on = t > 0 & t < 0.2;
%         selectIdx = cell(size(a));
% 
%         for k = 1 : length(a)
% 
%             % use contrast temporal stimuli
%             stimsForSelection = contains(a{k}.trials.events.trial_name, stimNm);
% 
%             % use HRF task, if not available, use other trials from spatialpattern
% 
%         %     stimsForSelection = contains(a{k}.trials.events.trial_name, {'pattern', 'HRFPATTERN', 'HRF', 'hrf'});
%         %     if max(stimsForSelection) == 0
%         %         % this patient has no HRF; use other trials from spatialpattern
%         %         stimsForSelection = contains(a{k}.trials.events.trial_name, {'SPARSITY', 'PLAID', 'CIRCULAR'});
%         %         stim_on = t > 0 & t < 0.2;
%         %     end
% 
%             mbb = mean(bb{k}(:,:,stimsForSelection),3);
%             llim = (mbb - (std(bb{k}(:,:,stimsForSelection),0,3)));
%             ulim = (mbb + (std(bb{k}(:,:,stimsForSelection),0,3)));
%             mbb_ci = cat(3, llim, ulim);
%             nEl = size(mbb,1); 
% 
%             % plot
%             figureName = sprintf('viselec_%s_%s_all', a{k}.sub, a{k}.ses);
%             figure('Name', figureName); nSubPlot = ceil(sqrt(nEl)); 
%             for el = 1:nEl
%                 subplot(nSubPlot,nSubPlot,el); hold on
%                 plotTitle = sprintf('%s W:%s B:%s ', a{k}.trials.channels.name{el}, a{k}.trials.channels.wangatlas{el}, a{k}.trials.channels.bensonatlas{el});        
%                 ecog_plotSingleTimeCourse(t, mbb(el,:), squeeze(mbb_ci(el,:,:)), [], plotTitle)       
%                 set(gcf, 'Position', [150 100 1500 1250]);
%             end
%             if saveFigure
%                 %printnice([], 0, saveLoc, figureName)
%                 saveas(gcf, fullfile(saveLoc, 'electrode_selection', figureName), 'png'); close;
%             end
% 
%             % FIRST SELECTION CRITERION:
%             select_idx1 = max(mbb(:,stim_on),[],2) > selTresh;
% 
%             % SECOND SELECTION CRITERION:
%             select_idx2 = zeros(size(select_idx1));
%             for el = 1:nEl
%                 if mean(mbb(el, stim_on)) > 0
%                     select_idx2(el) = 1;
%                 else
%                     select_idx2(el) = 0;
%                 end
%             end
% 
%             % EXCLUDE CHANNELS LABELED AS BAD
%             if isfield(summary(a{k}.trials.channels), 'status')
%                 select_idx3 = contains(a{k}.trials.channels.status, 'good');
%             else
%                 select_idx3 = ones(size(select_idx1));
%             end
% 
%             % EXCLUDE DEPTH ELECTRODES
%             select_idx4 = contains(a{k}.trials.channels.type, {'ecog', 'ECOG'});
% 
%             % combine criteria
%             select_idx = select_idx1 + select_idx2 + select_idx3 + select_idx4;
%             select_idx = (select_idx == 4);
% 
%             % plot selection
%             figureName = sprintf('viselec_%s_%s_selected', a{k}.sub, a{k}.ses);
%             figure('Name', figureName); nSubPlot = ceil(sqrt(nEl)); 
%             for el = 1:nEl
%                 if select_idx(el)
%                     subplot(nSubPlot,nSubPlot,el); hold on
%                     plotTitle = sprintf('%s W:%s B:%s ', a{k}.trials.channels.name{el}, a{k}.trials.channels.wangatlas{el}, a{k}.trials.channels.bensonatlas{el});
%                     ecog_plotSingleTimeCourse(t, mbb(el,:), squeeze(mbb_ci(el,:,:)), [], plotTitle)    
%                     set(gcf, 'Position', [150 100 1500 1250]);
%                 end
%             end
%             if saveFigure
%                 %printnice([], 0, saveLoc, figureName)
%                 saveas(gcf, fullfile(saveLoc, 'electrode_selection', figureName), 'png'); close;
%             end    
%             selectIdx{k} = select_idx;
%         end
 
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
    
end