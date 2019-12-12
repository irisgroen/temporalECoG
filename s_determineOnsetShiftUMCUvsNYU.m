% s_determineOnsetShiftUMCUvsNYU

% Purpose of this script
%
% Visual comparison of the broadband responses from V1 electrodes in the
% UMCU patient (sub-chaam) and NYU patients suggested a systematic offset in
% the onset of the visual response between the two sites. 
%
% After running several checks we concluded this must be due to a delay in
% the stimulus presentation relative to the trigger at UMCU. We decide to
% shift the stimuls onsets for the UMCU data by a fixed amount to bring the
% two datasets in alignment. 
%
% The shift is determined through cross-correlation of the evoked responses
% (not broadband) between sub-chaam and sub-....
%
% We perform the following steps:
% ...
%
%


% Load the data from patients with V1 coverage
subjectList = {'chaam', 'som692', 'som708'}; 
inputDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual/derivatives/ECoG';



