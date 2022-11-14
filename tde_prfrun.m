%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% DEMO script: extract and preprocess the ECoG PRF data and fit pRF model. %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 1: Load the ECoG PRF data and stimulus description

% Load and epoch the data
recomputeFlag = false;
tasks         = {'prf'};
epochTime     = [0 0.85];
saveStr       = 'prfdata';
sampleRate    = 512;
[data] = tde_getData(recomputeFlag, [], [], tasks, epochTime, sampleRate, [], saveStr);

% Compute the PRF timecourses
doPlots = false;
[data] = tde_computePRFtimecourses(data, [], [], doPlots);

% Load the stimulus apertures
stimName = fullfile(tdeRootPath, 'prf_apertures', 'bar_apertures.mat');
load(stimName, 'bar_apertures');

%% 2: Model fitting

% Fit the PRF time courses with analyzePRF
tr              = 1;
opt.hrf         = 1;
opt.maxpolydeg  = 0;
opt.xvalmode    = 0; 
opt.forcebounds = 1;
opt.display     = 'off';

doPlots = true;
[results] = tde_fitPRFs(data, bar_apertures, opt, doPlots);

 
