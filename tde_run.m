%% 1: Load the ECoG data and stimulus description

% load or (re)compute the processed data
reComputeFlag = false; 
[fulldata] = tde_getData(reComputeFlag);

% select epochs and channels, average trials within stimulus condition 
options = [];
options.average_elecs             = true;
options.elec_exclude_depth        = true;
options.doplots                   = false;
options.elec_selection_method     = 'splithalf';
options.areanames                 = 'V1';
%options.stimnames                 = {'CRF-1','CRF-2', 'CRF-3','CRF-4','CRF-5'}; %{'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4','ONEPULSE-5','ONEPULSE-6'}; %{'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};
[data, channels, stimnames, t, srate, options] = tde_selectData(fulldata, options);

% plot average response per stimulus for selected data
savePlot = 0; 
saveStr = [];%e.g. 'CRF';
tde_plotData(data, channels, t, options, savePlot, saveStr);

% generate stimulus timecourses
[stim_ts, stim_info] = tde_generateStimulusTimecourses(options.stimnames,t);

%% 2. Model fitting

% Define model(s)
modelfuns = tde_modelTypes();
modelfun = modelfuns([1 6 7 8 9 10 11 12 13 14]); 

% Define options
options.xvalmode = 0;      % 0 = none, 1 = stimulus leave-one-out
options.display  = 'off';  % 'iter' 'final' 'off'

LOADFITS = 1; % instead of fitting, load an existing saved model fit

if LOADFITS  
	% Load model fit(s)
    [params, pred] = tde_loadModelFits(modelfun, options);
else    
    % Compute model fit(s)
    [params, pred] = tde_doModelFits(modelfun, stim_ts, data, srate, options);
end

%% 3. Model evaluation

% Compute R2 and derived parameters
[results] = tde_evaluateModelFit(data, modelfun, params, pred, stim_info);

%% 4. Plot timecourses and fits

% Provide a directory to save figures (optional)
%saveDir = fullfile(analysisRootPath, 'figures', 'modelfits');
saveDir = [];
for ii = 1: length(results)
    tde_plotDataAndFits(results(ii), data, channels, stim_ts, stim_info, t, [], saveDir)
end

%tde_plotResiduals(results, data, channels, stim_ts, stim_info, t, [], saveDir)

%% 5. Plot derived and fitted parameters

saveDir = [];

% model parameters
%saveDir = fullfile(analysisRootPath, 'figures', 'modelparams');
tde_plotParams(results, channels, saveDir);close;

% model predictions (derived)
%saveDir = fullfile(analysisRootPath, 'figures', 'modelpredictions');
tde_plotDerivedPredictions(results,channels,2,1, saveDir);

%% UNDER DEVELOPMENT

% data params 
close all;
tde_plotDerivedParamsData(data,channels,t,stim_info)


%tde_computeDerivedParamsData(data,channels,t,stim_info);

