function [params, pred] = tde_doModelFits(modelfun, stim, data, channels, srate,t,options, saveDir, saveStr)

% Wrapper around tde_fitModel to run multiple models in a loop 

% <saveDir> path to save parameters and fits; if empty, results are not
%   saved (default)
% <saveName> string to add to the save filename, if results are saved
%   (default empty)

% Save options
if ~exist('saveStr', 'var'), saveName = []; end
if ~exist('saveDir', 'var'), saveDir = fullfile(analysisRootPath, 'results'); end

% Some formatting
if ~iscell(modelfun), modelfun = {modelfun}; end 

% Fit model(s)
params = []; pred = [];
if options.average_elecs
    preprocName = 'electrodeaverages';
else
    preprocName = 'individualelecs';
end

for ii = 1:size(modelfun,2)
    
    objFunction = modelfun{ii};
    [params{ii}, pred{ii}] = tde_fitModel(objFunction, stim, data, srate, options);
    
    % SAVE RESULTS
    if ~isempty(saveDir)

        if ~exist(saveDir, 'dir'); mkdir(saveDir); end
        if isempty(saveName)
            saveName = sprintf('%s_xvalmode%d_%s', func2str(objFunction), options.xvalmode, dataName);
        else
            saveName = sprintf('%s_xvalmode%d_%s_%s', func2str(objFunction), options.xvalmode, dataName, saveStr);
        end
        saveName = fullfile(saveDir, saveName);
        fprintf('[%s] Saving results to %s \n', mfilename, saveName);

        if exist(sprintf('%s.mat',saveName),'file')
            warning('[%s] Results file already exists! Writing new file with date-time stamp.',mfilename);
            saveName = sprintf('%s_%s', saveName, datestr(now,30));
            fprintf('[%s] Saving results to %s \n', mfilename, saveName);
        end
        save(saveName, 'pred', 'params', 'stim', 'data', 'channels','srate','t', 'options', 'objFunction');  
    end
end


