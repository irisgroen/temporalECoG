% script 

% load (1) or compute (0)
[epochs, channels, events] = tde_getData(0);

% next step, tde_selectElectrodes
% input: epochs, channels, events, taskstouse/trialstouse, inclusion criteria 

% next step, tde_fitModel
% input: model: function handle, file?
% output: model fits

% separate script?: tde_analyzePRF, ecog_analyzePRF?
% preferably NOT another save step, do on the fly? or add to channel table

% next step, tde_analyzeModelFits
% input: model fits (can be multiple?)
% outputs: plots, stats

