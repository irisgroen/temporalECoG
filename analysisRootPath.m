function rootPath = analysisRootPath()

% To reproduce figures reported in: Groen et al, (2022) Temporal dynamics of
% neural responses in human visual cortex (preprint:
% https://www.biorxiv.org/content/10.1101/2021.08.08.455547v3), update this
% path to point to the derivatives/Groen2022TemporalDynamics/ inside the
% Visual ECoG dataset, available on OpenNeuro:
% https://openneuro.org/datasets/ds004126.
%
% To run or reproduce the data analysis and generate new results, leave as
% is, and results and figures will be written to an 'analysis' folder
% that will be created inside the code repository (tdeRootPath).
%
% Iris Groen 2022

rootPath = fullfile(tdeRootPath, 'analysis');
if ~exist(rootPath, 'dir'), mkdir(rootPath); end


end

