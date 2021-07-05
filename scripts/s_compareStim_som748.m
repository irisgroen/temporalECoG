clear all;

subjects = {'som748'};
compute = false;

% get BAIR data
tasks = {'temporalpattern'};
saveStr = 'bair';
[data1] = tde_getData(compute, subjects, [], tasks, [], [], [], saveStr);

% select data
options = [];
options.stimnames = {'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
                     'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};
options.doplots = false;
[epochs1, channels1, t1, srate1, options1] = tde_selectData(data1, options);

% get SCENETEMPORAL data
tasks = {'sixcatloctemporal'};
saveStr = 'sixcat';
[data2] = tde_getData(compute, subjects, [], tasks, [], [], [], saveStr);

% select data
options = [];
options.stimnames = {'ONEPULSE-1','ONEPULSE-2', 'ONEPULSE-3','ONEPULSE-4', 'ONEPULSE-5','ONEPULSE-6',...
                     'TWOPULSE-1','TWOPULSE-2', 'TWOPULSE-3','TWOPULSE-4', 'TWOPULSE-5','TWOPULSE-6'};
options.doplots = false;
[epochs2, channels2, t2, srate2, options2] = tde_selectData(data2, options);

% create stim_info
[~,stim_info] = tde_generateStimulusTimecourses(options.stimnames,t1);
stim_info = stim_info(6:end,:);

%% compute recovery
w = 0.5;
metric = 'sum';
[ISIrecover1, ts1] = tde_computeISIrecovery(epochs1,t1,stim_info,srate1,w,[],metric);
[ISIrecover2, ts2] = tde_computeISIrecovery(epochs2,t2,stim_info,srate2,w,[],metric);

% keep only matching channels

[C,I1,I2] = intersect(channels1.name, channels2.name);
ISIrecover1 = ISIrecover1(:,I1);
ISIrecover2 = ISIrecover2(:,I2);
channels1 = channels1(I1,:);
channels2 = channels2(I2,:);

%% average within visual areas
areaNames = {'V123', 'V3ab', 'LOTO'};
areaNames = [];

[~, channels_m1, group_prob1] = groupElecsByVisualArea(channels1, 'probabilisticresample', areaNames);   
[m1, se1] = averageWithinArea(ISIrecover1, group_prob1, [], 1000);

[~, channels_m2, group_prob2] = groupElecsByVisualArea(channels2, 'probabilisticresample',  areaNames);   
[m2, se2] = averageWithinArea(ISIrecover2, group_prob2, [], 1000);

%% plot
figure('Position', [126         312        1206         458]);hold on;
channels = channels_m2;

conditionsOfInterest = {'ONEPULSE-4', 'TWOPULSE'};
stim_idx  = contains(stim_info.name, conditionsOfInterest);

cmap      = brewermap(height(channels)+2, 'RdBu');
cmap      = cmap(1:height(channels),:);
%cmap = [1 0 0; 0 0 1; 0 1 0];

x = stim_info.ISI(stim_idx)*1000; % in ms

chan_ind = [9 8 7 6 5 3 ];%height(channels):-1:1;
%chan_ind = [1 2 3];

subplot(1,2,1); cla; hold on

% Plot linear prediction 
h0 = line([x(1) x(end)], [1 1], 'LineStyle', ':', 'LineWidth', 2, 'color', [0 0 0]);
set(gca, 'xtick', x, 'xticklabelrotation', 45);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

for ii = chan_ind
    tde_plotPoints(m1(:,ii), se1(:,ii,:), x, 'ci', 0, [],[],cmap(ii,:));
end
ylim([0 2]);
xlim([-20 x(end)+20]);

% make legend
l = [];
for ii = 1:length(chan_ind)
    l{ii} = sprintf('%s (n=%d)', channels_m1.name{chan_ind(ii)}, channels_m1.number_of_elecs{chan_ind(ii)});
end
legend(l, 'location', 'southeast');

legend boxoff
xlabel('Stimulus interval (ms)');
ylabel('Recovery ratio');

title('BAIR (grayscale stimuli)');


subplot(1,2,2); cla; hold on

% Plot linear prediction 
h0 = line([x(1) x(end)], [1 1], 'LineStyle', ':', 'LineWidth', 2, 'color', [0 0 0]);
set(gca, 'xtick', x, 'xticklabelrotation', 45);
set(get(get(h0,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

for ii = chan_ind
    tde_plotPoints(m2(:,ii), se2(:,ii,:), x, 'ci', 0, [],[],cmap(ii,:));
end
ylim([0 2]);
xlim([-20 x(end)+20]);

% make legend
l = [];
for ii = 1:length(chan_ind)
    l{ii} = sprintf('%s (n=%d)', channels_m2.name{chan_ind(ii)}, channels_m2.number_of_elecs{chan_ind(ii)});
end
legend(l, 'location', 'southeast');

legend boxoff
xlabel('Stimulus interval (ms)');
ylabel('Recovery ratio');

title('SCENETEMPORAL (naturalistic images)');
set(findall(gcf,'-property','FontSize'),'FontSize',20)

