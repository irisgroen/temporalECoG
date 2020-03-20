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
%opts.areanames                 = 'V1';
%opts.stimnames                 = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5'};
%opts.stimnames                 = {'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6'};
%opts.stimnames                 = {'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};
[data, channels, stimnames, t, srate, opts] = tde_selectData(fulldata, opts);

% plot average response per stimulus for selected data
savePlot = 0;
tde_plotData(data, channels, t, opts, savePlot);

% generate stimulus timecourses
[stim_ts, stim_info] = tde_generateStimulusTimecourses(opts.stimnames,t);

%% 2. Model fitting

% define model(s)3
modelfuns = tde_modelTypes();
modelfun = modelfuns([1 2]); 

options          = [];
options.xvalmode = 0;      % 0 = none, 1 = stimulus leave-one-out
options.display  = 'off';  % 'iter' 'final' 'off

% define saveDir (optional)
saveDir  = fullfile(analysisRootPath, 'results');

% Fit model(s)
params = []; pred = [];
for ii = 1:size(modelfun,2)      
    [params{ii}, pred{ii}] = tde_fitModel(modelfun{ii}, stim_ts, data, srate, options, saveDir);
end

%% Or, load saved fits from disk

% define model(s)
modelfuns = tde_modelTypes();
modelfun = modelfuns([1]); 

% define model fitting options
options          = [];
options.xvalmode = 0;      % 0 = none, 1 = stimulus leave-one-out

% load in fits
params = []; pred = [];
for ii = 1:size(modelfun,2)
    if opts.average_elecs
        name = 'electrodeaverages';
    else
        name = 'individualelecs';
    end
	a = load(fullfile(saveDir, sprintf('%s_results_xvalmode%d_%s.mat', func2str(modelfun{ii}), options.xvalmode, name)));
    params{ii} = a.params;
    pred{ii} = a.pred;
end

%% 3. Model evaluation

% Compute R2 and derived parameters
group_inx = [1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3];
%group_inx = [1 1 1 1 1];
[results] = tde_evaluateModelFit(data, modelfun, params, pred, group_inx);

%% 4. Plot timecourses and fits

% Provide a directory so save figures (optional)
saveDir = fullfile(analysisRootPath, 'figures', 'modelfits');
saveDir = [];
tde_plotDataAndFits(results, data, channels, stim_ts, stim_info, t, [], saveDir)


saveDir = fullfile(analysisRootPath, 'figures', 'modelparams');
tde_plotParams(results, channels, saveDir);


saveDir = fullfile(analysisRootPath, 'figures', 'modelpredictions');
tde_plotDerivedPredictions(results,channels,2,1, saveDir);


