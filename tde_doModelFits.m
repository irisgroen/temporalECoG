function tde_doModelFits(modelfun,stim,data,channels,srate,t,stim_info,options,saveDir,saveStr)

% Wrapper around tde_fitModel to run multiple models in a loop. 

% <saveDir> path to save parameters and fits; if empty, results are not
%   saved (default)
% <saveStr> string to add to the save filename, if results are saved
%   (default empty)
%
% 2022 Iris Groen

% Save options
if ~exist('saveStr', 'var'), saveStr = []; end
if ~exist('saveDir', 'var') || isempty(saveDir), saveDir = fullfile(analysisRootPath, 'results'); end
if ~isfield(options,'fitaverage') || isempty(options.fitaverage), options.fitaverage = false; end

% Some formatting
if ~iscell(modelfun), modelfun = {modelfun}; end 

% Fit model(s)
if options.average_elecs
    preprocName = 'electrodeaverages';
else
    preprocName = 'individualelecs';
end

for ii = 1:size(modelfun,2)
    
    % FIT MODEL
    objFunction = modelfun{ii};
    
    if options.fitaverage
        if ~isfield(options,'nfits') || isempty(options.nfits)
            options.nfits = 10;
        end
        fprintf('[%s] Warning: options.fitaverage is %d; fitting each average %d times \n', mfilename, options.fitaverage, options.nfits);
        
        % Initialize
        params = [];
        pred = [];
        data_concat = [];
        channels_concat = [];
         
        % Fit multiple times to average
        for jj = 1:options.nfits
            fprintf('[%s] Starting fit on average %d of %d \n', mfilename, jj, options.nfits);
           
            [~, channels_av, group_prob] = groupElecsByVisualArea(channels, 'probabilisticresample', options.areanames);
            fun = @mean;
            numboot = 1; % no bootstrapping, just one assignment
            [data_av, ~, ~, n_elecs_selected] = averageWithinArea(data, group_prob, fun, numboot);
            [params_av, pred_av] = tde_fitModel(objFunction, stim, data_av, srate, options);
            
            % overwrite n_elecs_selected in channels table
            channels_av.number_of_elecs = n_elecs_selected;
            % add column to reflect number of fit
            channels_av.nfit = ones(height(channels_av),1) * jj;
            % remove subjects_name column because it doesn't match the selected elecs 
            channels_av = removevars(channels_av, {'subject_name'});
            % concatenate 
            params = cat(2,params,params_av);
            pred = cat(3,pred,pred_av);
            data_concat = cat(3,data_concat,data_av);
            channels_concat = vertcat(channels_concat, channels_av);
        end
        data = data_concat;
        channels = channels_concat;
    else
        % Fit once 
        [params, pred] = tde_fitModel(objFunction, stim, data, srate, options);
    end
    
    % SAVE RESULTS
    if ~isempty(saveDir)
        
        if ~exist(saveDir, 'dir'); mkdir(saveDir); end
        if isempty(saveStr)
            saveName = sprintf('%s_xvalmode%d_%s', func2str(objFunction), options.xvalmode, preprocName);
        else
            saveName = sprintf('%s_xvalmode%d_%s_%s', func2str(objFunction), options.xvalmode, preprocName, saveStr);
        end
        saveName = fullfile(saveDir, saveName);
        fprintf('[%s] Saving results to %s \n', mfilename, saveName);

        if exist(sprintf('%s.mat',saveName),'file')
            warning('[%s] Results file already exists! Writing new file with date-time stamp.',mfilename);
            saveName = sprintf('%s_%s', saveName, datestr(now,30));
            fprintf('[%s] Saving results to %s \n', mfilename, saveName);
        end
        save(saveName, 'stim', 'data', 'pred', 'params', 'channels','srate','t','stim_info','options', 'objFunction');  
    end
end


