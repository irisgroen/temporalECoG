function [params, pred] = tde_loadModelFits(modelfun, opts, resultsDir, saveStr)

% Loads saved model fits for multiple models in a loop
%
% 2020 Iris Groen

if ~exist('saveStr', 'var'), saveStr = []; end

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
    fname = fullfile(resultsDir, sprintf('%s_xvalmode%d_%s', modelName, opts.xvalmode, dataName));
    if ~isempty(saveStr)
        fname = sprintf('%s_%s', fname, saveStr);
    end     
	a = load(sprintf('%s.mat',fname));
    params{ii} = a.params;
    pred{ii} = a.pred;
end

end