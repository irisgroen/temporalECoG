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
opts.normalize_data  = false;
opts.average_elecs   = false;
[data2fit, channels, stimnames, t] = tde_selectData(data, [], opts);

% generate stimulus timecourses
[stim_ts] = tde_generateStimulusTimecourses(stimnames,t);

toc

%% fitting

ele = 54; 
smallData = data2fit(:,:,ele);

opts = [];
opts.srate = channels.sampling_frequency(ele);

%params:    [t1,   w, t2,   n,   sigma, shift, gain]
opts.x0   = [0.03, 0, 0.07, 1.5, 0.15, 0.06, 1];
opts.lb   = [0, 0, 0, 0, 0, 0, 0];
opts.ub   = [1, 1, 2, 10, 1, 1, inf];

opts.plb   = [.1, 0,.1, 1, .01, .01, .5];
opts.pub   = [.9, .5, 1, 3, .5, .5, 100];

tic
[results, pred] = tde_fitModel(@DNmodel, smallData, stim_ts, opts);
toc

figure;
subplot(2,1,1);plot(t,smallData, 'LineWidth', 3); title('data')
subplot(2,1,2);plot(t,pred, 'LineWidth', 3); title('tdefitmodel')


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



