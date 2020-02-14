%function tde_computePRFs

recomputeFlag = true;
subjects      = []; % will default to all subjects in subjectList.tsv
sessions      = []; % will default to all sessions per subject
tasks         = {'prf'};
description   = []; % will default to broadband
epochTime     = [-0.3 0.85];
sampleRate    = []; % will default to 512
saveStr       = 'prfdata';

[data] = tde_getData(recomputeFlag, subjects, sessions, tasks, description, epochTime, sampleRate, saveStr);


%end