%% 1: Load the ECoG data and stimulus description

% load or (re)compute the processed data
reComputeFlag = false; 
[fulldata] = tde_getData(reComputeFlag);

% select epochs and channels, average trials within stimulus condition 
opts = [];
opts.average_elecs             = true;
opts.elec_exclude_depth        = true;
opts.doplots                   = false;
opts.elec_selection_method     = 'splithalf';
opts.areanames                 = 'V1';
%opts.stimnames                 = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5'};
%opts.stimnames                 = {'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6'};
%opts.stimnames                 = {'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};
[data, channels, stimnames, t, srate, opts] = tde_selectData(fulldata, opts);

% plot average response per stimulus for selected data
savePlot = 0; 
saveStr = [];%'CRF';
tde_plotData(data, channels, t, opts, savePlot, saveStr);

% generate stimulus timecourses
[stim_ts, stim_info] = tde_generateStimulusTimecourses(opts.stimnames,t);

%% 2. Model fitting

% Define model(s)
modelfuns = tde_modelTypes();
modelfun = modelfuns([1 2]); 

% Define options
options          = [];
options.xvalmode = 0;      % 0 = none, 1 = stimulus leave-one-out
options.display  = 'off';  % 'iter' 'final' 'off

% Define saveDir (optional)
saveDir  = fullfile(analysisRootPath, 'results');

% Fit model(s)
params = []; pred = [];
for ii = 1:size(modelfun,2)      
    [params{ii}, pred{ii}] = tde_fitModel(modelfun{ii}, stim_ts, data, srate, options, saveDir);
end

% %% Or, load saved fits from disk
% 
% % Define model(s)
% modelfuns = tde_modelTypes();
% modelfun  = modelfuns([1 2]); 
% xvalmode  = 0;
% 
% [params, pred] = tde_loadSavedModelFit(modelfun, xvalmode, opts);

%% 3. Model evaluation

% Compute R2 and derived parameters
[results] = tde_evaluateModelFit(data, modelfun, params, pred, stim_info);

%% 4. Plot timecourses and fits

% Provide a directory so save figures (optional)
saveDir = fullfile(analysisRootPath, 'figures', 'modelfits');
tde_plotDataAndFits(results, data, channels, stim_ts, stim_info, t, [], saveDir)

%% 5. Plot derived and fitted parameters

% model parameters
saveDir = fullfile(analysisRootPath, 'figures', 'modelparams');
tde_plotParams(results, channels, saveDir);

% model predictions (derived)
saveDir = fullfile(analysisRootPath, 'figures', 'modelpredictions');
tde_plotDerivedPredictions(results,channels,2,0, saveDir);


%% UNDER DEVELOPMENT

% data params 
%saveDir = fullfile(analysisRootPath, 'figures', 'data');
%close all;
shift = params{1}(6,:);
tde_computeDerivedParamsData(data,channels,t,shift,stim_info);

