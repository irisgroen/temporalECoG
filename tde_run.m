%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% DEMO script: extract and preprocess the ECoG data and fit temporal models.%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1: Load the ECoG data 

% Load or (re)compute the processed data
reCompute = false; 
[data_full] = tde_getData(reCompute);

% Select epochs and channels, average trials within stimulus condition
options = [];
options.doplots = false;
[data_selection, channels_selection, t, srate, options] = tde_selectData(data_full, options);

% Generate stimulus timecourses
[stim_ts, stim_info] = tde_generateStimulusTimecourses(options.stimnames,t);

% Sort and average electrodes
options.average_elecs = true;
options.normalize_data = false;  % boolean
options.areanames = {'V1','V2','V3','V3a','V3b','LO1','LO2','TO','IPS'}; % exclude hV4, combine TO1 and TO2
[data, channels, se] = tde_prepareData(data_selection, channels_selection, options);

% Plot average response per stimulus for selected data
tde_plotData(data, channels, t, options,1);

%% 2. Model fitting

% Define model(s)
modelfuns = tde_modelTypes();
modelfun = modelfuns([1]);

% Define options
options.xvalmode = 1;      % 0 = none, 1 = stimulus leave-one-out
options.display  = 'off';  % 'iter' 'final' 'off'
options.algorithm = 'bads';
options.fitaverage = false;
options.nfits = 1000; % if fit average 

% Compute model fit(s); data and fits will be saved to 'analysis/results' folder
tde_doModelFits(modelfun, stim_ts, data, channels, srate, t, stim_info, options);

%% 3. Model evaluation

% Load data and fits
modelfun = @TTCSTIG19;
xvalmode = 0;
datatype = 'electrodeaverages';
[D] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% Compute R2 and derived parameters
objFunction = modelfun;
[results] = tde_evaluateModelFit(D,1);

%% 4. Plot timecourses and fits

% Provide a directory to save figures (optional)
saveDir = fullfile(analysisRootPath, 'figures');

% Plot multiple model predictions (superimposed)
tde_plotDataAndFits(results, D.data, D.channels, D.stim, D.stim_info, D.t, saveDir, {'ONEPULSE', 'TWOPULSE', 'CRF'})
%tde_plotResiduals(results, D.data, D.channels, D.stim, D.stim_info, D.t, saveDir, {'ONEPULSE', 'TWOPULSE', 'CRF'})

%% 5. Plot derived and fitted parameters

% model parameters
saveDir = [];%fullfile(analysisRootPath, 'figures', 'modelparams');
tde_plotParams(results, D.channels, saveDir);%close;

% model predictions (derived)
saveDir = fullfile(analysisRootPath, 'figures', 'modelpredictions');
tde_plotDerivedPredictions(results,D.channels,1,1,saveDir);
tde_plotDerivedPredictions(results,D.channels,1,0,saveDir);
tde_plotDerivedPredictions(results,D.channels,2,1,saveDir);
tde_plotDerivedPredictions(results,D.channels,2,0,saveDir);

%% UNDER DEVELOPMENT

% data params 
%close all;
tde_plotDerivedParamsData(pred{1},channels,t,stim_info, [],0);
tde_plotDerivedParamsData(D.data,D.channels,D.t,D.stim_info, {'V1', 'V2', 'V3'},0);
tde_plotDerivedParamsModel(params{1},modelfun{1}, channels,t,{'V1'});


