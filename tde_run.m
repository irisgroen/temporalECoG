%% 1: Load the ECoG data 

% Load or (re)compute the processed data
reComputeFlag = false; 
[data_full] = tde_getData(reComputeFlag);

% Select epochs and channels, average trials within stimulus condition
options = [];
options.doplots = false;
[data_selection, channels_selection, t, srate, options] = tde_selectData(data_full, options);

% Generate stimulus timecourses
[stim_ts, stim_info] = tde_generateStimulusTimecourses(options.stimnames,t);

% Sort and average electrodes
options.average_elecs = false;
options.normalize_data = false;  % boolean
[data, channels, se] = tde_prepareData(data_selection, channels_selection, options);

% Plot average response per stimulus for selected data
tde_plotData(data, channels, t, options);

%% 2. Model fitting

% Define model(s)
modelfuns = tde_modelTypes();
modelfun = modelfuns([10]); 

% Define options
options.xvalmode = 0;      % 0 = none, 1 = stimulus leave-one-out
options.display  = 'off';  % 'iter' 'final' 'off'
options.algorithm = 'bads';

LOADFITS = 0; % instead of fitting, load an existing saved model fit
saveStr = [];%'sixROIs';

if LOADFITS  
	% Load model fit(s)
    [params, pred] = tde_loadModelFits(modelfun, options, [], saveStr);
else    
    % Compute model fit(s)
    [params, pred] = tde_doModelFits(modelfun, stim_ts, data, channels, srate, t, stim_info, options);
end

%% 3. Model evaluation

% Compute R2 and derived parameters
objFunction = modelfun;
[results] = tde_evaluateModelFit(data, objFunction, params, pred, stim_info);

%% 4. Plot timecourses and fits

% Provide a directory to save figures (optional)
saveDir = fullfile(analysisRootPath, 'figures');

% % Plot data and prediction for each model individually
% for ii = 1 : length(results)
%     tde_plotDataAndFits(results(ii), data2fit, channels2fit, stim_ts, stim_info, t, saveDir,  {'ONEPULSE', 'TWOPULSE', 'CRF'});    %tde_plotResiduals(results(ii), data, channels, stim_ts, stim_info, t, saveDir)
% end


% Plot multiple model predictions (superimposed)
tde_plotDataAndFits(results, data, channels, stim_ts, stim_info, t, saveDir, {'ONEPULSE', 'TWOPULSE', 'CRF'})

%% HACK to plot areas superimposed
tmp = results(1);
for ii = 5%2:9
    tmp.pred(:,:,1) = data2fit(:,:,ii);
    tde_plotDataAndFits(tmp, data, channels, stim_ts, stim_info, t, [],{'ONEPULSE'})
end
set(gcf, 'Position', [121         622        1800         363]);

tmp = results(1);
for ii = 5%2:9
    tmp.pred(:,:,1) = tmp.pred(:,:,ii);
    tde_plotDataAndFits(tmp, results(1).pred, channels2fit, stim_ts, stim_info, t, [], {'ONEPULSE'})
end
set(gcf, 'Position', [121         622        1800         363]);

%% 5. Plot derived and fitted parameters

%saveDir = [];

% model parameters
saveDir = fullfile(analysisRootPath, 'figures', 'modelparams');
tde_plotParams(results, channels, []);%close;

% model predictions (derived)
saveDir = fullfile(analysisRootPath, 'figures', 'modelpredictions');
tde_plotDerivedPredictions(results,channels,2,1, saveDir);

%% UNDER DEVELOPMENT

% data params 
%close all;
tde_plotDerivedParamsData(data,channels,t,stim_info, [],0);

%tde_plotDerivedParamsData(pred{5},channels,t,stim_info, {'V1', 'V2', 'V3'}, 0)
%tde_computeDerivedParamsData(data,channels,t,stim_info);

tde_plotDerivedParamsModel(params{1},modelfun{1}, channels,t,{'V1'});


