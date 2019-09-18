% script 

% load (1) or compute (0)
tic
[data] = tde_getData(1);
toc
% next step, tde_selectElectrodes
% input: epochs, channels, events, taskstouse/trialstouse, inclusion criteria 

% next step, tde_generateStimulusTimecourses
% input: events
% output: time series per event

% next step, tde_fitModel
% input: model: function handle, data, eventCodes, stimulusTimeCourses, opts
% output: model fits

% -- which models?
% ----- DN (flavors: uniphasic, biphasic, fixed exponent or not)
% ----- DN cascade?
% ----- Two temporal channels (flavors: HH, Stigliani 1, Stigliani 2)
% ----- DN-like models:
% --------- Heeger 1993

% separate script?: tde_analyzePRF, ecog_analyzePRF?
% preferably NOT another save step, do on the fly? or add to channel table

% next step, tde_analyzeModelFits
% input: model fits (can be multiple?)
% outputs: plots, stats

