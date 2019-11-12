% script 

tic

% load or (re)compute the processed data
reComputeFlag = false; 
[data] = tde_getData(reComputeFlag);

% select epochs and channels, average trials within stimulus condition 
data_opts = [];
data_opts.doplots         = false;
data_opts.normalize_data  = false;
data_opts.average_elecs   = false;
[data2fit, channels, stimnames, t, srate] = tde_selectData(data, [], data_opts);

% generate stimulus timecourses
[stim_ts] = tde_generateStimulusTimecourses(stimnames,t);

toc

%% fitting

% define electrode (temporary)
ele = 54; 
smallData = data2fit(:,:,ele);

% define model
modelType = 'TTC'; 

modelfun = str2func(sprintf('%smodel', modelType));

tic
[results, pred] = tde_fitModel(modelfun, smallData, stim_ts, srate);
toc

figure;
subplot(2,1,1);plot(t,smallData, 'LineWidth', 3); title('data')
subplot(2,1,2);plot(t,pred, 'LineWidth', 3); title('tdefitmodel')

% -- which models?
% ----- DN (flavors: uniphasic, biphasic, fixed exponent or not)
% ----- DN cascade?
% ----- Two temporal channels (flavors: HH, Stigliani 1, Stigliani 2)
% ----- DN-like models:
% --------- Heeger 1993

% next step, tde_analyzeModelFits
% input: model fits 
% ---> put deriveParams and R2 in here (rather than in tde_modelFit)
% outputs: plots, stats

% function tde_plotModelFits 


% How about the PRFs --> separate pipeline (like this one), read in fits
% from file (e.g. add to channel table), prf_getData, prf_selectData,
% prf_fitModel



