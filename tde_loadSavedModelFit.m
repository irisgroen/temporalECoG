function [params, pred] = tde_loadSavedModelFit(modelfun, xvalmode, opts)

if ~iscell(modelfun), modelfun = {modelfun}; end 

% define saveDir (optional)
saveDir  = fullfile(analysisRootPath, 'results');

% load in fits
params = []; pred = [];
for ii = 1:size(modelfun,2)
    if opts.average_elecs
        name = 'electrodeaverages';
    else
        name = 'individualelecs';
    end
	a = load(fullfile(saveDir, sprintf('%s_xvalmode%d_%s.mat', func2str(modelfun{ii}), xvalmode, name)));
    params{ii} = a.params;
    pred{ii} = a.pred;
end

end