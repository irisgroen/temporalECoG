function rootPath = analysisRootPath()

% Users of the temporalECoG repository should update this path to point to
% where they want the processed data and figures to be written to.

%rootPath = '/Volumes/server/Projects/BAIR/Papers/TemporalDynamicsECoG';
rootPath = fullfile(tdeRootPath, 'analysis');
if ~exist(rootPath, 'dir'), mkdir(rootPath); end


end

