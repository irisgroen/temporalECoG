%% 1: Load the ECoG PRF data and stimulus description

% Load and epoch the data
recomputeFlag = true;
tasks         = {'prf'};
epochTime     = [-0.2 0.6];
saveStr       = 'prfdata';
sampleRate    = 512;

[data] = tde_getData(recomputeFlag, [], [], tasks, epochTime, sampleRate, [], saveStr);

% Compute the PRF timecourses
doPlots = true;

[data] = tde_computePRFtimecourses(data, [], [], doPlots);

% Load the stimulus apertures
stimName = fullfile(tdeRootPath, 'prf_apertures', 'bar_apertures.mat');
load(stimName, 'bar_apertures');

%% 2: Model fitting

% Fit the PRF time courses with analyzePRF
tr             = 1;
opt.hrf        = 1;
opt.maxpolydeg = 0;
opt.xvalmode   = 0; 
opt.display    = 'off';

doPlots = true;
[results] = tde_fitPRFs(data, bar_apertures, opt, doPlots);

%% 3: Model evaluation? E.g. Summarize R2 within participants?


%% 4. Plot timecourse fits and estimated prfs - individual elecs
% to facilitate comparison with tde model fits


% Add some summary plots comparing e.g. benson and analyzePRF?
% Separately for each subject + across subjects