function [stim] = tde_generateStimulusTimecourses(stimNames,t, stimInfo)
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

% %% old jing version of stim ts generation
% stim  = zeros(nStim,length(t));
% contrasts = [0.0405, 0.0902, 0.2105, 0.3203, 1.0000];
% durs      = [0.016667, 0.033333, 0.066667, 0.13333, 0.26667, 0.53333];
%
% % CONTRAST STIMULI --------------------------------------------------------
% stim(1 : 5, t > 0 & t<=0.5) = 1;
% for k = 1 : 5, stim(k, :) = stim(k, :) .* contrasts(k); end
% 
% % INCREASING DURATIONS ----------------------------------------------------
% for k = 1 : 6, stim(k + 5, (t>0) & (t <= durs(k))) = 1; end
% 
% % INCREASING ISI ----------------------------------------------------------
% stim(12 : nStim, t > 0 & t <= durs(4)) = 1;
% for k = 1 : 6
%     t_start = durs(4) + durs(k);
%     t_end   = durs(4) * 2 + durs(k);
%     stim(11 + k, t > t_start & t <= t_end) = 1;
% end
% 
% stim = stim'; % put time dimension first
% figure;
% for ii = 1:length(stimNames)
%     subplot(nStim,1,ii); plot(t,stim(:,ii),'LineWidth', 2); ylim([0 1.2]);
% end   


end