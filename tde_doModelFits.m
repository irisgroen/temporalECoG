function [params, pred] = tde_doModelFits(modelfun, stim_ts, data, srate, options, saveDir)

% Wrapper around tde_fitModel to run multiple models in a loop 

if ~exist('saveDir', 'var') || isempty(saveDir) 
    saveDir = fullfile(analysisRootPath, 'results'); 
end

if ~iscell(modelfun), modelfun = {modelfun}; end 

% Fit model(s)
params = []; pred = [];
if options.average_elecs
    saveStr = 'electrodeaverages';
else
    saveStr = 'individualelecs';
end

for ii = 1:size(modelfun,2)      
    [params{ii}, pred{ii}] = tde_fitModel(modelfun{ii}, stim_ts, data, srate, options, saveDir, saveStr);
end
