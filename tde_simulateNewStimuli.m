function [stim, stim_info] = tde_simulateNewStimuli(t,nStim)

% Generate stimulus time courses to use for generating new model
% predictions after having fitted a model to the experimental data.
%
% 2022 Iris Groen

% CRF
% generate stimuli with nStim contrast values
stimCRF = zeros(length(t),nStim);
nameCRF = cell(nStim,1);
durationCRF = zeros(nStim,1)+0.5;
isiCRF = zeros(nStim,1);
contrastCRF = (1/nStim:1/nStim:nStim/nStim)';

stim_on = t>0 & t<=0.5;
for ii = 1:nStim
    stimCRF(stim_on,ii) = ii/nStim;
    nameCRF{ii} = sprintf('CRF-sim%d', ii);
end

% DUR
% generate stimuli with nStim durations between 0 and 0.533;
stimDUR = zeros(length(t),nStim);
nameDUR = cell(nStim,1);
durationDUR = zeros(nStim,1);
isiDUR = zeros(nStim,1);
contrastDUR = ones(nStim,1);

for ii = 1:nStim
    t_off = ii/nStim * 0.533;
    stim_on = t>0 & t<=t_off;
    stimDUR(stim_on,ii) = 1;
    durationDUR(ii) = round(t_off,3);
    nameDUR{ii} = sprintf('ONEPULSE-sim%d', ii);
end

% ISI
% generate stimuli with nStim ISIs between 0 and 0.533;
stimISI = zeros(length(t),nStim);
nameISI = cell(nStim,1);
durationISI = ones(nStim,1)*0.133;
isiISI = zeros(nStim,1);
contrastISI = ones(nStim,1);
% generate first pulse
pulse1_toff = 0.133;
pulse1_on = t>0 & t<=pulse1_toff;
stimISI(pulse1_on,:) = 1;
ISIs = linspace(0,0.533,nStim);
% generate second pulse
for ii = 1:nStim
    this_ISI = ISIs(ii);%ii/nStim * 0.533;
    pulse2_ton = pulse1_toff + this_ISI;
    pulse2_toff = pulse2_ton + 0.133;
    pulse2_on = t>pulse2_ton & t<=pulse2_toff;
    stimISI(pulse2_on,ii) = 1;
    isiISI(ii) = round(this_ISI,3);
    nameISI{ii} = sprintf('TWOPULSE-sim%d', ii);
end

% concatenate conditions
stim = horzcat(stimCRF, stimDUR, stimISI);
name = vertcat(nameCRF, nameDUR, nameISI);
duration = vertcat(durationCRF, durationDUR, durationISI);
ISI =  vertcat(isiCRF, isiDUR, isiISI);
contrast = vertcat(contrastCRF, contrastDUR, contrastISI);

% generate new stim_info table
stim_info = table(name, duration, ISI, contrast);    

end