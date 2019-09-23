% script 

% note: sub-beilen not included at the moment because of inconsistencies in
% the formatting of the events.tsv files, need to sync with Gio.

% load (1) or compute (0)
tic
[data] = tde_getData(0);
toc

% select epochs and channels 
tic
[data2] = tde_selectData(data,1);
toc

% tde_averageEpochs/concatData
[data3] = tde_combineData(data2,0);



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

