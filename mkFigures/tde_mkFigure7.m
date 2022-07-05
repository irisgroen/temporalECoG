% Figure 7 : 
% Delayed normalization explains multiple temporal dynamics in multiple visual areas
%
% 2022 Iris Groen

modelfuns = {@LINEAR,@LINEAR_RECTF,@LINEAR_RECTF_EXP,@LINEAR_RECTF_EXP_NORM,@LINEAR_RECTF_EXP_NORM_DELAY};

xvalmode = 1;
datatype = 'individualelecs';

for ii = 1:length(modelfuns)
    [d(ii)] = tde_loadDataForFigure(modelfuns{ii}, xvalmode, datatype);
end

[results] = tde_evaluateModelFit(d,0);
numboot = 10000;

nModels = length(modelfuns);
cmap = flipud(brewermap(nModels,'Spectral'));

%% Panel A (here figure 1)

% Subplot positions: % [left bottom width height]
posa = [0.05 0.10 0.95 0.85];

figure(1); clf; hold on
set(gcf, 'position',  get(0, 'screensize'));

% Cross-validated R2 for model build up in each area
modelind = [1 2 3 4 5];
R2 = nan(length(modelind),height(d(1).channels));

for ii = 1:size(R2,1)
    R2(ii,:) = results(modelind(ii)).R2.concat_all;
end

[~, channels, group_prob] = groupElecsByVisualArea(d(1).channels, 'probabilisticresample');   
[m, se] = averageWithinArea(R2, group_prob, [], numboot);

subplot('position', posa);  cla; hold on
x = 1:height(channels);
tde_plotBars(m,se,x,cmap);

l = {'Linear', '+ Rectification', '+ Exponentiation', '+ Normalization', '+ Delay', '+ Delay (fixed IRF)'};
legend(l, 'location', 'northeast');
legend('boxoff')
set(gca, 'xtick', 1:height(channels), 'xticklabel', channels.name)
set(gca, 'ylim', [0 1], 'xlim', [0 size(m,2)+1]);
box off
ylabel('Cross-validated R^{2}', 'Interpreter','tex'); 

set(findall(gcf,'-property','FontSize'),'FontSize',24)

%% Panel B (here figure 2)
plotSE = 1;

% Subplot positions: % [left bottom width height]
posa = [0.04 0.58 0.4 0.4];
posb = [0.5 0.58 0.4 0.4];
posc = [0.04 0.08 0.4 0.4];
posd = [0.5 0.08 0.4 0.4];

figure(2); clf; hold on
set(gcf, 'position',  get(0, 'screensize'));

group_prob = group_prob(:,1); % V1 electrodes

D = d(1);
t = D.t;
stim_ts = D.stim;
stim_info = D.stim_info;
nModels = length(modelfuns); % omit DN

% 1. temporal summation

% Find stimulus index
stim_idx = find(contains(stim_info.name, 'ONEPULSE'));
x = stim_info.duration(stim_idx)*1000; % in ms

% Generate new stimuli based on params
nStimNew = 50;
[stim2, stim_info2] = tde_simulateNewStimuli(t,nStimNew);
stim_info2.duration = stim_info2.duration*1000; % convert to ms
stim_idx2 = find(contains(stim_info2.name, 'ONEPULSE'));
x2 = stim_info2.duration(stim_idx2); % in ms

% Predict model responses for new stimuli
pred2 = nan([nModels size(stim2(:,stim_idx2)) height(D.channels)]);
srate = D.channels.sampling_frequency(1);
for kk = 1:nModels
    DD = d(kk);
    for ii = 1:size(DD.params,2)
        prm = DD.params(:,ii);
        [~, pred2(kk,:,:,ii)] = DD.objFunction(prm, [], stim2(:,stim_idx2), srate);      
    end
end

% Determine time index over which to compute summary statistics
t_idx = t>0 & t<1;

% Concatenate data and prediction (in order to make sure probabilistic
% assignment of electrodes to areas is done the same way for both).
dat = D.data(t_idx,stim_idx,:);
for kk = 1:nModels
    dat = cat(2,dat,squeeze(pred2(kk,t_idx, :,:)));
end

% Compute sum across stim_on window
m_conc = squeeze(sum(dat,1)); 
[m_conc, se_conc] = averageWithinArea(m_conc, group_prob, @median, numboot);

% Plot
subplot('position', posa);cla; hold on

% Plot prediction
nStim = length(stim_idx);
M = m_conc(nStim+1:end);
M = reshape(M, [nStimNew nModels]);
SE = se_conc(nStim+1:end,:);
SE = reshape(SE, [nStimNew nModels 2]);
for kk = 1:nModels
    m = M(:,kk);
    se = SE(:,kk,:);
    if plotSE == 1
        tde_plotPoints(m, se, x2, 'ci', 1,[],50, cmap(kk,:));
    else
        tde_plotPoints(m, [], x2, 'ci', 1,[],50, cmap(kk,:));
    end
end

% Plot data
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'errbar', 1,[],50);

l = {'Linear', '+ Rectification', '+ Exponentiation', '+ Normalization', '+ Delay'};

% Format axes
ticklabelsX = num2str(x); ticklabelsX(2:5,:) = ' ';
set(gca, 'xtick', x, 'xticklabel', ticklabelsX);
xlabel('Stimulus duration (ms)'); ylabel({'Summed broadband', 'power (0-1s)'}); 
lgd = legend([l {'Neural data'}], 'location', 'none');
lgd.Position = [0.45 0.65 0.01 0.1];
%posa = [0.04 0.58 0.4 0.4];

%title ('Temporal summation');

legend('boxoff')
axis tight
axis square
xlim([0 550]);
ylim([0 1.3]);
set(findall(gcf,'-property','FontSize'),'FontSize',24)

% 2. Repetition suppression

% Find stimulus index
stim_idx = find(contains(stim_info.name, {'ONEPULSE-5', 'TWOPULSE'}));
x = stim_info.ISI(stim_idx)*1000; 

% Compute recovery per electrode
srate = D.channels.sampling_frequency(1);
[m,ts] = tde_computeISIrecovery(D.data,D.t,D.stim_info,srate, [], [], 'max');

% Generate new TWOPULSE stimuli based on params
nStimNew = 50;
[stim2, stim_info2] = tde_simulateNewStimuli(t,nStimNew);
stim_idx2 = find(contains(stim_info2.name, 'TWOPULSE'));
stim2 = stim2(:,stim_idx2);
stim_info2 = stim_info2(stim_idx2,:);
x2 = stim_info2.ISI*1000; 

% Add the ONEPULSE-4 condition
stim_idx1 = find(contains(stim_info.name, 'ONEPULSE-4'));
stim2 = cat(2,D.stim(:, stim_idx1),stim2);
stim_info2 = [stim_info(stim_idx1,[1 3:5]); stim_info2];

% Predict model responses for new stimuli
pred2 = nan([nModels size(stim2) height(D.channels)]);
srate = D.channels.sampling_frequency(1);
for kk = 1:nModels
    DD = d(kk);
    for ii = 1:size(DD.params,2)
        prm = DD.params(:,ii);
        [~, pred2(kk,:,:,ii)] = DD.objFunction(prm, [], stim2, srate); 
    end
end

% Compute recovery per electrode
m2 = []; ts2 = [];
for kk = 1:nModels
    [m2(kk,:,:), ts2(kk,:,:,:)] = tde_computeISIrecovery(squeeze(pred2(kk,:,:,:)),t,stim_info2,srate);
end

% Concatenate data and prediction (in order to make sure probabilistic
% assignment of electrodes to areas is done the same way for both).
for kk = 1:nModels
    m = cat(1,m,squeeze(m2(kk,:,:)));
end

% Compute average parameter values within groups
[m_conc, se_conc] = averageWithinArea(m, group_prob, @median, numboot);

% Plot
subplot('position', posb); cla; hold on

% Plot prediction
nStim = length(stim_idx);
M = m_conc(nStim+1:end);
M = reshape(M, [nStimNew nModels]);
SE = se_conc(nStim+1:end,:);
SE = reshape(SE, [nStimNew nModels 2]);
for kk = 1:nModels
    m = M(:,kk);
    se = SE(:,kk,:);
    if plotSE == 1
        tde_plotPoints(m, se, x2, 'ci', 0,[],50, cmap(kk,:));
    else
        tde_plotPoints(m, [], x2, 'ci', 0,[],50, cmap(kk,:));
    end
end

% Plot data
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'errbar', 0,[],50);

% Format axes
ticklabelsX = num2str(x); ticklabelsX(2:6,:) = ' ';
set(gca, 'xtick', x, 'xticklabel', ticklabelsX);

ylim([0 1.2]);
xlim([-20 x(end)+20]);
xlabel('Stimulus interval (ms)'); ylabel('Ratio second / first stimulus'); 
axis tight
axis square
%title ('Repetition suppression');

% 3-4 Contrast dynamics

% Find stimulus index
stim_idx = find(contains(stim_info.name, 'CRF'));
x = stim_info.contrast(stim_idx) * 100; 

% Generate new stimuli based on params
nStimNew = 50;
[stim2, stim_info2] = tde_simulateNewStimuli(t,nStimNew);
stim_idx2 = find(contains(stim_info2.name, 'CRF'));
x2 = stim_info2.contrast(stim_idx2) * 100; 

% Predict model responses for new stimuli
pred2 = nan([nModels size(stim2(:,stim_idx2)) height(D.channels)]);
srate = D.channels.sampling_frequency(1);
for kk = 1:nModels
    DD = d(kk);
    for ii = 1:size(DD.params,2)
        prm = DD.params(:,ii);
        [~, pred2(kk,:,:,ii)] = DD.objFunction(prm, [], stim2(:,stim_idx2), srate);      
    end
end

% Determine time index over which to compute summary statistics
t_idx = t>0.05 & t<1.0;

% Concatenate data and prediction (in order to make sure probabilistic
% assignment of electrodes to areas is done the same way for both).
dat = D.data(t_idx,stim_idx,:);
for kk = 1:nModels
    dat = cat(2,dat,squeeze(pred2(kk,t_idx, :,:)));
end

% 3. Time to peak
subplot('position', posc); cla; hold on

% Compute peak across trial window
[~,I] = max(dat,[],1);
T = t(t_idx);
[m_conc, se_conc] = averageWithinArea(squeeze(T(I)), group_prob,  @median, numboot);

m_conc = m_conc * 1000; % convert to ms
se_conc = se_conc * 1000;


% Plot prediction
nStim = length(stim_idx);
M = m_conc(nStim+1:end);
M = reshape(M, [nStimNew nModels]);
SE = se_conc(nStim+1:end,:);
SE = reshape(SE, [nStimNew nModels 2]);
for kk = 1:nModels
    m = M(:,kk);
    se = SE(:,kk,:);
    if plotSE == 1
        tde_plotPoints(m, se, x2, 'ci', 0,[],50, cmap(kk,:));
    else
        tde_plotPoints(m, [], x2, 'ci', 0,[],50, cmap(kk,:));
    end
end

% Plot data
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'errbar', 0, [], 50);

% Format axes
l = get(gca, 'YLim'); ylim([50 l(2)]);
ticklabelsX = num2str(x); ticklabelsX(2:4,:) = ' ';
set(gca, 'xlim', [0 105], 'xtick', x, 'xticklabel', ticklabelsX);
xlabel('Contrast (%)'); ylabel('Time-to-peak (ms)'); 
axis tight
axis square
%title ('Contrast time-to-peak');

% 4. Ratio sustained/transient
subplot('position', posd); cla; hold on
T = D.t(t_idx);

% Smooth time courses to get better estimates of max and offset response levels
% Apply same amount of smoothing to both data and model predictions
dat2 = dat;
for ii = 1:size(dat,2)
    for jj = 1:size(dat,3)
        dat2(:,ii,jj) = smooth(dat(:,ii,jj),150);
    end
end

[M] = max(dat2,[],1); % value at peak
t_off = (T == 0.5); % offset timepoint
O = dat2(t_off,:,:); % value at offset
R = squeeze(O./M); % divide value at offset with value at peak

[m_conc, se_conc] = averageWithinArea(R, group_prob, @median, numboot);

% Plot prediction
nStim = length(stim_idx);
M = m_conc(nStim+1:end);
M = reshape(M, [nStimNew nModels]);
SE = se_conc(nStim+1:end,:);
SE = reshape(SE, [nStimNew nModels 2]);
for kk = 1:nModels
    m = M(:,kk);
    se = SE(:,kk,:);
    if plotSE == 1
        tde_plotPoints(m, se, x2, 'ci', 0,[],50, cmap(kk,:));
    else
        tde_plotPoints(m, [], x2, 'ci', 0,[],50, cmap(kk,:));
    end
end

% Plot data
nStim = length(stim_idx);
m = m_conc(1:nStim);
se = se_conc(1:nStim,:);
tde_plotPoints(m, se, x, 'errbar', 0, [], 50);

% Format axes
ticklabelsX = num2str(x); ticklabelsX(2:4,:) = ' ';
set(gca, 'xlim', [0 105], 'xtick', x, 'xticklabel', ticklabelsX);
set(gca, 'ylim', [0 1]);
xlabel('Contrast (%)'); ylabel('Proportion offset/peak'); 
%title ('Contrast ratio transient/sustained');
axis tight
axis square

set(findall(gcf,'-property','FontSize'),'FontSize',24)


%%%% obsolete %%%%

% %% Plot
% 
% % 1. Contrast response function
% subplot('position', posc); cla; hold on
% 
% Compute sum across stim_on window
% m_conc = squeeze(sum(dat,1)); 
% [m_conc, se_conc] = averageWithinArea(m_conc, group_prob, @median, numboot);

% % Plot data
% nStim = length(stim_idx);
% m = m_conc(1:nStim);
% se = se_conc(1:nStim,:);
% tde_plotPoints(m, se, x, 'errbar', 1, [], 50);
% 
% % Plot prediction
% M = m_conc(nStim+1:end);
% M = reshape(M, [nStimNew nModels]);
% SE = se_conc(nStim+1:end,:);
% SE = reshape(SE, [nStimNew nModels 2]);
% for kk = 1:nModels
%     m = M(:,kk);
%     se = SE(:,kk,:);
%     if plotSE == 1
%         tde_plotPoints(m, se, x2, 'ci', 1,[],50, cmap(kk,:));
%     else
%         tde_plotPoints(m, [], x2, 'ci', 1,[],50, cmap(kk,:));
%     end
% end
% 
% % Format axes
% ticklabelsX = num2str(x); ticklabelsX(2:4,:) = ' ';
% set(gca, 'xlim', [0 105], 'xtick', x, 'xticklabel', ticklabelsX);
% xlabel('Contrast (%)'); ylabel('Summed broadband power (0-1s)'); 
% axis tight
% axis square
% title ('Contrast response function');
% set(findall(gcf,'-property','FontSize'),'FontSize',24)
