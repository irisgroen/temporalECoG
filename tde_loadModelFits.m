function [params, pred] = tde_loadModelFits(modelfun, opts, resultsDir)

% Loads saved model fits for multiple models in a loop 

if ~exist('resultsDir', 'var') || isempty(resultsDir) 
    resultsDir = fullfile(analysisRootPath, 'results'); 
end

if ~iscell(modelfun), modelfun = {modelfun}; end 

% load in fits
if opts.average_elecs
    dataName = 'electrodeaverages';
else
    dataName = 'individualelecs';
end

params = []; pred = [];
for ii = 1:size(modelfun,2)
    modelName = func2str(modelfun{ii});
    fprintf('[%s] Loading fits for %s for model %s \n', mfilename, dataName, modelName);
	a = load(fullfile(resultsDir, sprintf('%s_xvalmode%d_%s.mat', modelName, opts.xvalmode, dataName)));
    params{ii} = a.params;
    pred{ii} = a.pred;
end

end