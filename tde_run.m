%% 1: Load the ECoG data and stimulus description

% load or (re)compute the processed data
reComputeFlag = false; 
[data] = tde_getData(reComputeFlag);

% select epochs and channels, average trials within stimulus condition 
opts = [];
opts.doplots         = false;
opts.average_trials  = true;
opts.normalize_data  = false;
opts.average_elecs   = true;
opts.sort_channels   = true;
[data2fit, channels, stimnames, t, srate] = tde_selectData(data, [], opts);

% generate stimulus timecourses
[stim_ts, stim_info] = tde_generateStimulusTimecourses(stimnames,t);

% plot average response per stimulus for selected data
tde_plotData(data2fit, channels, t, opts);

%% 2. Model fitting

% define subset of data (temporary)
tmpdata = data2fit;%(:,:,1);

% define model(s)
modelfuns = tde_modelTypes();
modelfun = modelfuns([1 3]); 

% define model fitting options
options = struct();
options.xvalmode = 1;      % 0 = none, 1 = stimulus leave-one-out
options.display  = 'iter'; % 'iter' 'final' 'off

% define saveDir (optional)
saveDir = '/Volumes/server/Projects/BAIR/Papers/TemporalDynamicsECoG/results';

% fit model(s)
for ii = 1:size(modelfun,2)
    tic
    [params{ii}, pred{ii}] = tde_fitModel(modelfun{ii}, stim_ts, tmpdata, srate, options, saveDir);
    toc
end
        
%% 3. Model fit evaluation

% compute R2 and derived parameters
[results] = tde_evaluateModelFit(tmpdata, modelfun, params, pred);

% plot timecourses and fits
tde_plotDataAndFits(results, tmpdata, channels, stim_ts, stim_info, t)

% plot r2, derived params and fitted params
tde_plotFittedAndDerivedParams(results, channels);

%%
% -- which models?
% ----- DN (flavors: uniphasic, biphasic, fixed exponent or not)
% ----- DN cascade?
% ----- Two temporal channels (flavors: HH, Stigliani 1, Stigliani 2)
% ----- DN-like models:
% --------- Heeger 1993

% How about the PRFs --> separate pipeline (like this one), read in fits
% from file (e.g. add to channel table), prf_getData, prf_selectData,
% prf_fitModel



