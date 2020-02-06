%% 1: Load the ECoG data and stimulus description

% load or (re)compute the processed data
reComputeFlag = false; 
[data] = tde_getData(reComputeFlag);

% select epochs and channels, average trials within stimulus condition 
opts = [];
opts.doplots              = true;
opts.average_elecs        = true;
opts.elec_max_thresh      = 0.5;
opts.elec_exclude_depth   = true;
[data2fit, channels, stimnames, t, srate] = tde_selectData(data, [], opts);

% plot average response per stimulus for selected data
tde_plotData(data2fit, channels, t, opts);

% generate stimulus timecourses
[stim_ts, stim_info] = tde_generateStimulusTimecourses(stimnames,t);

%% 2. Model fitting

% define model(s)
modelfuns = tde_modelTypes();
modelfun = modelfuns([1 3 4 5 6 7]); 
%modelfun = modelfuns([1 7]); 

% define model fitting options
options          = [];
options.xvalmode = 1;      % 0 = none, 1 = stimulus leave-one-out
options.display  = 'off';  % 'iter' 'final' 'off

% define saveDir (optional)
saveDir  = fullfile(analysisRootPath, 'results');
LOADFITS = 1;

% Fit or load model(s)
params = []; pred = [];

for ii = 1:size(modelfun,2)
    
    if LOADFITS
        if opts.average_elecs
            name = 'electrodeaverages';
        else
            name = 'individualelecs';
        end
        a = load(fullfile(saveDir, sprintf('%s_results_xvalmode%d_%s.mat', func2str(modelfun{ii}), options.xvalmode, name)));
        params{ii} = a.params;
        pred{ii} = a.pred;
    
    else        
        tic
        [params{ii}, pred{ii}] = tde_fitModel(modelfun{ii}, stim_ts, data2fit, srate, options, saveDir);
        toc
    end
end

%% 3. Model evaluation

% Compute R2 and derived parameters
group_inx = [1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 3];
[results] = tde_evaluateModelFit(data2fit, modelfun, params, pred, group_inx);

%% 4. Plot timecourses and fits

% Provide a directory so save figures (optional)
saveDir = fullfile(analysisRootPath, 'figures', 'modelfits');

tde_plotDataAndFits(results, data2fit, channels, stim_ts, stim_info, t, [], saveDir)


saveDir = fullfile(analysisRootPath, 'figures', 'modelparams');

tde_plotParams(results, channels, saveDir);


saveDir = fullfile(analysisRootPath, 'figures', 'modelpredictions');

tde_plotDerivedPredictions(results, channels,1,0, saveDir);


