% script 

% load (1) or compute (0)
tic
[data] = tde_getData(0);
toc

% select electrodes, select epochs (or split up in separate scripts?)
stimNames = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5',...
'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};

[data] = tde_selectData(data, stimNames, baselineTime);


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

