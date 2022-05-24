function [stim, stimInfo] = tde_generateStimulusTimecourses(stimNames,t, stimInfo)

% Generates stimulus timecourses based on stimulus name and epoch duration,
% provided in a stimInfo table; will load stimInfo table in the code
% repository if no stimInfo input argument is provided.
%
% 2022 Iris Groen

% <subjectList>
if ~exist('stimInfo', 'var') || isempty(stimInfo)
    fprintf('[%s] Loading stimulus info \n',mfilename);
    stimInfo_fname = fullfile(tdeRootPath, 'stiminfo.tsv');
    stimInfo = readtable(stimInfo_fname, 'FileType', 'text');
end

fprintf('[%s] Generating stimulus timecourses \n',mfilename);

condition = stimInfo.name;
duration  = stimInfo.duration;
ISI       = stimInfo.ISI;
contrast  = stimInfo.contrast;

nStim = length(stimNames);
stim  = zeros(length(t),nStim);

for ii = 1:length(stimNames)
    inx = contains(condition, stimNames{ii});
    stim(t>0 & t<=duration(inx),ii) = contrast(inx);
    if ISI(inx) > 0
        t_pulse1off = duration(inx)+ISI(inx);
        stim(t>t_pulse1off & t<=(t_pulse1off+duration(inx)),ii) = contrast(inx);
    end
end
fprintf('[%s] Done! \n',mfilename);

end