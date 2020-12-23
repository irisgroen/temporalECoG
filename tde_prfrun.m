%% 1: Load the ECoG PRF data and stimulus description

% Load and epoch the data
recomputeFlag = false;
tasks         = {'prf'};
epochTime     = [-0.2 0.6];
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

%% 3: Model evaluation? E.g. Summarize R2 within participants?

% Add some summary plots comparing e.g. benson and analyzePRF?
% Separately for each subject + across subjects
% Compare with Ken's numbers


% %Following is the script to load the pRF data (without decimation). I
% think it's easier than exploring in the finder. "gaussianmode" can be
% 'og' (One Gaussian), 'dog' (Difference of Gaussians), or 'gs' (now we
% call Large-Field Suppression).
% 
% subjectList = {'som726'};
%  
% average        ='runs';
% smoothingMode  ='none';
% smoothingN     = [];
% prfmodel       ='linear';
% gaussianmode   ='gs';
%  
% opts = [];
% opts.average        = average;
% opts.smoothingMode  = smoothingMode;
% opts.smoothingN     = smoothingN;
% opts.prfmodel       = prfmodel;
% opts.gaussianmode   = gaussianmode;
% opts.issave         = false;
% opts.compute        = false; 
% 
% opts.targetBAND     ='bbS';
% prf_params_bb = ecog_prf_analyzePRF(subjectList, opts);
 
