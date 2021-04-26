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

dataName = fullfile(dataPath, sprintf('%s.mat', dataName));

if exist(dataName, 'file')
    data = load(dataName);
else
    error('Could not locate datafile')
end

if ~isfield(data.options, 'fitaverage')
    data.options.fitaverage = 0;
end
end