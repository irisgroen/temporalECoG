function [data] = tde_loadDataForFigure(modelfun, xvalmode, dataType, dataPath, dataStr)

if ~exist('dataStr', 'var')
    dataStr = [];
end

if ~exist('dataPath', 'var') || isempty(dataPath)
    dataPath = fullfile(analysisRootPath, 'results');
end

if isempty(dataStr)
    dataName = sprintf('%s_xvalmode%d_%s', func2str(modelfun), xvalmode, dataType);
else
    dataName = sprintf('%s_xvalmode%d_%s_%s', func2str(modelfun), xvalmode, dataType, dataStr);
end

data = load(fullfile(dataPath, dataName));
   
end