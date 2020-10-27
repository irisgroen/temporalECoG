%% 1: Load the ECoG data 

% Load or (re)compute the processed data
reComputeFlag = false; 
[fulldata] = tde_getData(reComputeFlag);

% Select epochs and channels, average trials within stimulus condition 
options.doplots = false;
[data, channels, t, srate, options] = tde_selectData(fulldata, options);

% Generate stimulus timecourses
[stim_ts, stim_info] = tde_generateStimulusTimecourses(options.stimnames,t);

% Sort and average electrodes
options.average_elecs = false;
options.area_names = []; %{'V1', 'V2', 'V3', 'V3ab', 'LOTO', 'IPS'};
[data2fit, channels2fit, se] = tde_prepareData(data, channels, options);

% Plot average response per stimulus for selected data
tde_plotData(data2fit, channels2fit, t, options);

%% 2. Model fitting

% Define model(s)
modelfuns = tde_modelTypes();
modelfun = modelfuns([6 10]); 

% Define options
options.xvalmode = 0;      % 0 = none, 1 = stimulus leave-one-out
options.display  = 'off';  % 'iter' 'final' 'off'

LOADFITS = 1; % instead of fitting, load an existing saved model fit
saveStr = [];%'sixROIs';

if LOADFITS  
	% Load model fit(s)
    [params, pred] = tde_loadModelFits(modelfun, options, [], saveStr);
else    
    % Compute model fit(s)
    [params, pred] = tde_doModelFits(modelfun, stim_ts, epochs, srate, options);
end

%% 3. Model evaluation

% Compute R2 and derived parameters
[results] = tde_evaluateModelFit(epochs, modelfun, params, pred, stim_info);

%% 4. Plot timecourses and fits

% Provide a directory to save figures (optional)
saveDir = fullfile(analysisRootPath, 'figures', 'modelfits_sixROIs');
% Plot data and prediction for each model individually
for ii = 1 : length(results)
    tde_plotDataAndFits(results(ii), epochs, channels2fit, stim_ts, stim_info, t, saveDir)
    %tde_plotResiduals(results(ii), data, channels, stim_ts, stim_info, t, saveDir)
end

% Plot multiple model predictions (superimposed)
tde_plotDataAndFits(results([2 1]), epochs, channels2fit, stim_ts, stim_info, t, saveDir)

%% 5. Plot derived and fitted parameters

%saveDir = [];

% model parameters
saveDir = fullfile(analysisRootPath, 'figures', 'modelparams');
tde_plotParams(results, channels2fit, saveDir);%close;

% model predictions (derived)
saveDir = fullfile(analysisRootPath, 'figures', 'modelpredictions');
tde_plotDerivedPredictions(results,channels2fit,2,1, saveDir);

%% UNDER DEVELOPMENT

% data params 
%close all;
tde_plotDerivedParamsData(data2fit,channels2fit,t,stim_info, {'V1', 'V2', 'V3'}, 0);

%tde_plotDerivedParamsData(pred{5},channels,t,stim_info, {'V1', 'V2', 'V3'}, 0)
%tde_computeDerivedParamsData(data,channels,t,stim_info);

tde_plotDerivedParamsModel(params{1},modelfun{1}, channels2fit,t,{'V1'});
