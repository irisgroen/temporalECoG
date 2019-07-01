
% Check whether we have the ECoG_utils repository on the path
if ~exist('ecog_matchChannels', 'file')
    tbUse ECoG_utils;
end

addpath(genpath('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal'));
addpath(genpath('/Volumes/server/Projects/Temporal_integration/DN_2018_data/code'))

stimNm = {'CRF-1','CRF-2', 'CRF-3','CRF-4', 'CRF-5',...
'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};

saveFigure = 0;
saveLoc = fullfile(dn_ctrst_RootPath, 'dataFigures');
    
%% load the data

fileName = 'preproc_epoched_visualelecs_notshifted.mat';
load(fullfile('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess', fileName));

%% 
t = a{1}.trials.time;

% extract evoked activity
ev = {};
for k = 1 : length(a)
    ev{k} = a{k}.trials.evoked;
end

% baseline correction

% use the firsrt 200ms before stimulus onset as the baseline
base_range = (t >= -0.2 & t < 0);

for k = 1 : length(a)
    %m_base{k} = squeeze(median(mean(ev{k}(:, base_range, :), 2), 3));
    %ev{k}     = ev{k}./m_base{k};
    baseline = mean(ev{k}(:,base_range,:),2);
    ev{k} = ev{k} - baseline;
end

%%

%% single trials

k = 1; %648
test = squeeze(ev{k}(ecog_matchChannels('MO01', a{k}.trials),:,find(contains(a{k}.trials.events.trial_name, stimNm))));
figure;plot(t, test); title([a{k}.sub ' single trials']);

k = 5; %708
test = squeeze(ev{k}(ecog_matchChannels('O02', a{k}.trials),:,find(contains(a{k}.trials.events.trial_name, stimNm))));
figure;plot(t, test); title([a{k}.sub ' single trials']);

k = 10; %chaam
test = squeeze(ev{k}(ecog_matchChannels('Oc18', a{k}.trials),:,find(contains(a{k}.trials.events.trial_name, stimNm))));
figure;plot(t, test); title([a{k}.sub ' single trials']);

k = 10; %chaam
test = squeeze(ev{k}(ecog_matchChannels('Oc12', a{k}.trials),:,find(contains(a{k}.trials.events.trial_name, stimNm))));
figure;plot(t, test); title([a{k}.sub ' single trials']);

%% trial averages
figure; hold on
subs = [1 10];
elecs = {'MO01', 'Oc12'};
colors = {'r', 'b'};
lnames = [];
for ee = 1:length(elecs)
    k = subs(ee);
    el = ecog_matchChannels(elecs{ee}, a{k}.trials);
    st = find(contains(a{k}.trials.events.trial_name, stimNm));
    temp = squeeze((ev{k}(el,:,st)));
    timecourse = mean(temp,2);
    ci = [timecourse-std(temp,0,2) timecourse+std(temp,0,2)];
    %ci = [];
    ecog_plotSingleTimeCourse(t,timecourse,ci,colors{ee}, [], 'Evoked response');
    lnames{ee} = [a{k}.sub ' ' elecs{ee} ' ' a{k}.trials.channels.bensonatlas{el}];
end
legend(lnames);

%% trial averages

% 708 vs chaam
subs = [5 10];
elecs = {'O02', 'Oc18'}; %V1 
%elecs = {'O03', 'Oc12'}; %V2

% 648 vs chaam
subs = [1 10];
elecs = {'MO01', 'Oc12'};

% 692 vs chaam
subs = [4 10];
elecs = {'MIO03', 'Oc18'}; %V1 
elecs = {'BO01', 'Oc12'}; %V2


colors = {'r', 'b'};
lnames = [];
onsetIndices = [];
figure; hold on

for ee = 1:2
    k = subs(ee); disp(a{k}.sub);
    el = ecog_matchChannels(elecs{ee}, a{k}.trials);
    st = find(contains(a{k}.trials.events.trial_name, stimNm));
    temp = squeeze((ev{k}(el,:,st)));
    timecourse = mean(temp,2);
%         if ee == 1
%             timecourse = -timecourse;
%         end
    %ci = [timecourse-std(temp,0,2) timecourse+std(temp,0,2)];
   
    ci = [];
    ecog_plotSingleTimeCourse(t,timecourse,ci,colors{ee}, [], 'Evoked response');
    lnames{ee} = [a{k}.sub ' ' elecs{ee} ' ' a{k}.trials.channels.bensonatlas{el}];
    
    difftest = timecourse';
    %difftest = [0 diff(test)];
    %p1 = plot(t,difftest,[colors{ee} ':'], 'LineWidth', 2);
    thresh = 0.15*max(abs(difftest));
    ind = find((difftest > thresh | difftest < -thresh) & t>0);
    %p2 = plot(t(ind(1)), test(ind(1)), '.k','MarkerSize', 50, 'LineStyle', 'none');
    p3 = plot(t(ind(1)), difftest(ind(1)), '.k','MarkerSize', 50, 'LineStyle', 'none');
    %p1.Annotation.LegendInformation.IconDisplayStyle = 'off';
    %p2.Annotation.LegendInformation.IconDisplayStyle = 'off';
    p3.Annotation.LegendInformation.IconDisplayStyle = 'off';

    onsetIndices(ee) = ind(1);
    
end
shiftInSamples = diff(onsetIndices);
shiftInSeconds = shiftInSamples * 1/a{1}.trials.fsample;
disp([num2str(shiftInSamples) ' ' num2str(1000*shiftInSeconds)]);

ecog_plotSingleTimeCourse(t-shiftInSeconds,timecourse,ci,[colors{ee} ':'], [], 'Evoked response');
p4 = plot(t(ind(1)-shiftInSamples), difftest(ind(1)-shiftInSamples), 'ok','MarkerSize', 15, 'LineStyle', 'none');
p4.Annotation.LegendInformation.IconDisplayStyle = 'off';

lnames = [lnames {[lnames{ee} ' shifted']}];
legend(lnames);

