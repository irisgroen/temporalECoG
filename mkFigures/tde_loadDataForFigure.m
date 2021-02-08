function [data] = tde_loadDataForFigure(modelfun, xvalmode, dataType, dataPath)

if ~exist('dataPath', 'var') || isempty(dataPath)
    dataPath = fullfile(analysisRootPath, 'results');
end

dataName = sprintf('%s_xvalmode%d_%s', func2str(modelfun), xvalmode, dataType);
data = load(fullfile(dataPath, dataName));     
   
end