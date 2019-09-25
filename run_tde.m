% script 

% note: sub-beilen not included at the moment because of inconsistencies in
% the formatting of the events.tsv files, need to sync with Gio.
tic

% load (0) or (re)compute (1)
reComputeFlag = false; 
[data] = tde_getData(reComputeFlag);

% select epochs and channels, average trials within stimulus condition 
doPlotsFlag = false; 
[data2fit, channels, stimnames, t] = tde_selectData(data, doPlotsFlag);

% generate stimulus timecourses
[stim_ts] = tde_generateStimulusTimecourses(stimnames,t);

toc

%% 

opts = [];
opts.normalize_data  = true;
opts.average_elecs   = true;

[out] = tde_fitModel(@dn_DNmodel, data2fit, stim_ts, t, channels, opts);

% next step, tde_fitModel
% input: model: function handle, data, eventCodes, stimulusTimeCourses, opts
% output: model fits
% -- which models?
% ----- DN (flavors: uniphasic, biphasic, fixed exponent or not)
% ----- DN cascade?
% ----- Two temporal channels (flavors: HH, Stigliani 1, Stigliani 2)
% ----- DN-like models:
% --------- Heeger 1993

% next step, tde_analyzeModelFits
% input: model fits (can be multiple?)
% outputs: plots, stats

% function tde_plotModelFits 
% How about the PRFs --> separate pipeline (like this one), read in fits
% from file (e.g. add to channel table), prf_getData, prf_selectData,
% prf_fitModel



