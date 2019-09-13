
[epochs, channels, events] = tde_getData(recomputeflag 0/1);


% next script, tde_prepElectrodes
% input: epochs, channels, events, taskstouse/trialstouse, inclusion criteria 

% next script, tde_fitModel
% input: model: function handle, file?
% output: model fits

% next script, tde_analyzeModelFits
% input: model fits (can be multiple?)
% outputs: plots, stats

% separate script: tde_analyzePRF
% preferably NOT another save step, do on the fly? or add to channel table
