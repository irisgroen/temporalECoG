function [stim, stimInfo] = tde_generateStimulusTimecourses(stimNames,t, stimInfo)
% Generates a timecourse based on stimulus name and epoch duration TO DO
% make it actually parse the stimName, OR, even better construct afresh
% from duration and ISI columns in events table (and stimulus mat file for
% contrasts??)

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