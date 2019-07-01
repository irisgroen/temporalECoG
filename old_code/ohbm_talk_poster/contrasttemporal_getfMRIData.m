tbUse bairanalysis;


setenv('SUBJECTS_DIR', '/Volumes/server/Freesurfer_subjects/');

projectDir           = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
modelType            = 'threeTasksUpsampledSurface';
tasks                = ['spatialpattern','temporalpattern','spatialobject'];
conditionsOfInterest = ["crf", "onepulse", "twopulse"];

subjectList = {'wlsubj052', 'wlsubj053', 'wlsubj054'};

B = []; SE = [];
for ii = 1:length(subjectList)
    
    subject = subjectList{ii};

    [meanBeta,betasSE,GLMconditions] = bidsSummarizeGLMDenoisebyArea (projectDir , subject, modelType, [], tasks, conditionsOfInterest);
    
    trialIndex = contains(GLMconditions, {'CRF', 'ONEPULSE', 'TWOPULSE'});
    B(ii,:,:) = meanBeta(:,trialIndex);
    SE(ii,:,:) = betasSE(:,trialIndex);
end

%[meanBeta,betasSE , GLMconditions] = bidsSummarizeGLMDenoisebyArea (projectDir , subject, modelType,...
    %session,tasks, conditionsOfInterest, makeFigures, saveFigures);
    
saveName = 'preproc_BAIRfmri';
save(fullfile('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess', saveName), 'B','SE');
    
%mfmri_2fit = [];
%sfmri_2fit = [];