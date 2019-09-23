function [stim] = tde_generateStimulusTimecourses(stimNames,t)


% 
contrasts = [0.0625 0.125 0.25 0.5 1]; 
durs      = [0.016667, 0.033333, 0.066667, 0.13333, 0.26667, 0.53333];

nStim = length(stimNames);
stim  = zeros(nStim, length(t));

% CONTRAST STIMULI --------------------------------------------------------
stim(1 : 5, t > 0 & t<=0.5) = 1;
for k = 1 : 5, stim(k, :) = stim(k, :) .* contrasts(k); end

% INCREASING DURATIONS ----------------------------------------------------
for k = 1 : 6, stim(k + 5, (t>0) & (t <= durs(k))) = 1; end

% INCREASING ISI ----------------------------------------------------------
stim(12 : nStim, t > 0 & t <= durs(4)) = 1;
for k = 1 : 6
    t_start = durs(4) + durs(k);
    t_end   = durs(4) * 2 + durs(k);
    stim(11 + k, t > t_start & t <= t_end) = 1;
end

end