% script 

% note: sub-beilen not included at the moment because of inconsistencies in
% the formatting of the events.tsv files, need to sync with Gio.
tic

% load (0) or (re)compute (1)
reComputeFlag = false; 
[data] = tde_getData(reComputeFlag);

% select epochs and channels, average trials within stimulus condition 
opts = [];
opts.doplots         = false;
opts.normalize_data  = true;
opts.average_elecs   = false;
[data2fit, channels, stimnames, t] = tde_selectData(data, [], opts);

% generate stimulus timecourses
[stim_ts] = tde_generateStimulusTimecourses(stimnames,t);

srate = 1/median(diff(t)); % samples per second
toc

%% 
% To call fitting function, we need:
%   1. the objective function (model form and type of error)
%   2. data
%   3. stimuli
%   4. starting values and bounds for parameters (there should be defaults
%               for each type of model)
%   5. sample rate of the data

ele = 54; 
smallData = data2fit(:,:,ele);

opts = [];
opts.srate = srate;
% biphasic
% opts.x0   = [0.03, 0, 0.07, 1.5, 0.15, 0.06, 1];
% opts.lb   = [0, 0, 0, 0, 0, 0, 0];
% opts.ub   = [1, 1, 1, 10, 1, 1, 1];
% uniphasic
opts.x0   = [0.03, 0.07, 1.5, 0.15, 0.06, 1];
opts.lb   = [0, 0, 0, 0, 0, 0];
opts.ub   = [1, 1, 10, 1, 1, 1];

[results] = tde_fitModel(@DNmodel, smallData, stim_ts, opts);

[results] = tde_fitModel(@DNmodel, data2fit, stim_ts, t, channels, opts);

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



