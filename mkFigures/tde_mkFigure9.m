% tde_mkFigure 9

% Load data and fits
modelfun = @DN;%;LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datastr = 'fitaverage1000bads';
datatype = 'individualelecs';
D = tde_loadDataForFigure(modelfun, xvalmode, datatype, [], datastr);

% Compute derived parameters
[results] = tde_evaluateModelFit(D);

% Select electrodes and compute averages 
if D.options.fitaverage
    % this means electrodes were averaged within area and fitted multiple times
    [data, data_se, channels] = averageMultipleFits(D.data, D.channels, @mean);
    [pred, pred_se] = averageMultipleFits(D.pred, D.channels, @mean);
else
    % this means we have just one fit to all individual electrodes 
    [~, channels, group_prob] = groupElecsByVisualArea(D.channels, 'probabilisticresample');   
    [data, data_se] = averageWithinArea(D.data, group_prob, @mean, 1000);
    [pred, pred_se] = averageWithinArea(D.pred, group_prob, @mean, 1000);
end
t = D.t;
stim_ts = D.stim;
stim_info = D.stim_info;

%% Plot specs
% Subplot positions: % [left bottom width height]

posa1 = [0.05  0.75 0.9 0.22];
posa2 = [0.05  0.50 0.9 0.22];
posb =  [0.05  0.1 0.25 0.3];
posc =  [0.375 0.1 0.25 0.3];
posd =  [0.7   0.1 0.25 0.3];

figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

% Panel A: time courses comparison V1 and V3b

conditionsOfInterest = {'CRF'};
timepointsOfInterest = [-0.1 1.2];
areasOfInterest      = {'V1', 'V3b'};

stim_idx  = contains(stim_info.name, conditionsOfInterest);
t_idx     = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);
chan_idx1 = find(strcmp(channels.name, areasOfInterest{1}));
chan_idx2 = find(strcmp(channels.name, areasOfInterest{2}));

% Top: Data

subplot('position', posa1); hold on

s = stim_ts(t_idx,stim_idx);
c1 = data(t_idx,stim_idx,chan_idx1);
c2 = data(t_idx,stim_idx,chan_idx2);
c1 = c1/max(c1(:,end));
c2 = c2/max(c2(:,end));

hs = plot(s(:), 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(c1(:), 'k-', 'linewidth', 2);
hp = plot(c2(:), 'Color', ones(1,3)*0.5, 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xtick',1:size(c1,1):size(c1,2)*size(c1,1), 'xticklabel',[]);
set(gca, 'ytick', [0 1]);
box off,  axis tight

ylim([-0.2 1.5]);
xlim([-20 length(s(:)) + 20]);

ylabel(sprintf('Neural response \n (normalized)')); 
legend(areasOfInterest, 'location', 'northeast');
legend('boxoff');

% Bottom: Model

subplot('position', posa2); hold on

s = stim_ts(t_idx,stim_idx);
c1 = pred(t_idx,stim_idx,chan_idx1);
c2 = pred(t_idx,stim_idx,chan_idx2);
c1 = c1/max(c1(:,end));
c2 = c2/max(c2(:,end));

hs = plot(s(:), 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(c1(:), 'r-', 'linewidth', 2);
hp = plot(c2(:), 'Color', [1 0.5 0.5], 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xtick',1:size(c1,1):size(c1,2)*size(c1,1), 'xticklabel', stim_info.contrast(stim_idx)*100, 'xticklabelrotation', 45);
set(gca, 'ytick', [0 1]);
box off,  axis tight

ylim([-0.2 1.5]);
xlim([-20 length(s(:)) + 20]);

xlabel('Contrast (%)'); ylabel(sprintf('Model prediction \n (normalized)')); 
legend(areasOfInterest, 'location', 'northeast');
legend('boxoff');

set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% Panel B: c50 summarized across areas

% Data

% Determine time index over which to compute sum
t_idx = t>0.05 & t<1.0;

% Concatenate data and prediction (in order to make sure probabilistic
% assignment of electrodes to areas is done the same way for both).
d = D.data(t_idx,stim_idx,:);

%d_s = smooth(d(:),10);
%d = reshape(d_s, size(d));

% Compute sum across stim_on window
sumd = squeeze(max(d,[],1));
%sumd = squeeze(sum(d,1));

c50_data = nan(1,height(D.channels));
fprintf('[%s] Fitting Naka Rushton function to each dataset...\n',mfilename);
for ii = 1:height(D.channels)
    [~,c50_data(ii)] = fitNakaRushton(stim_info.contrast(stim_idx)*100,sumd(:,ii),0);
end

if D.options.fitaverage
    [m_data, se_data, channels] = averageMultipleFits(c50_data, D.channels, @mean);
else
    [m_data, se_data] = averageWithinArea(c50_data, group_prob, @mean, 1000);
end

% Model
c50_model = results.derived.params(4,:);

if D.options.fitaverage
    [m_model, se_model, channels] = averageMultipleFits(c50_model, D.channels, @mean);
else
    [m_model, se_model] = averageWithinArea(c50_model, group_prob, @mean, 1000);
end

subplot('position', posb); cla; hold on
x = 1:height(channels);
tde_plotPoints(m_model', se_model, x, 'ci', 0, [], 40, 'r');
tde_plotPoints(m_data', se_data, x, 'errbar', 0, [], 40, 'k');

set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
ylabel('C50 (% contrast)');
xlabel('Visual area');

xlim([0 max(x)+1]); ylim([0 100]);
legend('model', 'data', 'Location','NorthWest');legend boxoff
   
%% Panel C: scaling error across areas

% Data
p = D.pred(:,stim_idx,:);
cScale_data = nan(2,height(D.channels));
cScale_model = nan(2,height(D.channels));
fprintf('[%s] Regressing low and high contrast conditions for each dataset...\n',mfilename);
for ii = 1:height(D.channels)
    %X = d(:,1,ii)/max(d(:,1,ii)); % 6.25% contrast
    %y = d(:,5,ii)/max(d(:,5,ii)); % 100% contrast
    X = d(:,1,ii); % 6.25% contrast
    y = d(:,5,ii); % 100% contrast
    
    %mdl = fitlm(X,y,'Intercept',false);
    mdl = fitlm(X,y);
    cScale_data(1,ii) = sum(abs(mdl.Residuals.Raw)); 
    cScale_data(2,ii) = mdl.Rsquared.Ordinary;

   % tmp = pdist([X y]', 'correlation');
   % cScale_data(1,ii) = tmp; 
   
    %X = p(:,1,ii)/max(p(:,1,ii)); % 6.25% contrast
    %y = p(:,5,ii)/max(p(:,5,ii)); % 100% contrast
    X = p(:,1,ii); % 6.25% contrast
    y = p(:,5,ii); % 100% contrast
    
   % tmp = pdist([X y]', 'correlation');
   % cScale_model(1,ii) = tmp; 
	%mdl = fitlm(X,y,'Intercept',false);
	mdl = fitlm(X,y);

    cScale_model(1,ii) = sum(abs(mdl.Residuals.Raw)); 
    cScale_model(2,ii) = mdl.Rsquared.Ordinary;

end
%figure;plot(t(t_idx),X,t(t_idx),y);
%title(tmp);
%tmp = mdl.Residuals.Raw;
%figure;plot(t(t_idx),X,t(t_idx),y,t(t_idx),tmp)
%legend('low','high','resid')
%title(mdl.Rsquared.Ordinary);

if D.options.fitaverage
    [m_data, se_data, channels] = averageMultipleFits(cScale_data, D.channels, @mean);
else
    [m_data, se_data] = averageWithinArea(cScale_data, group_prob, @mean, 1000);
end

% Model
%cScale_model = results.derived.params([5 6],:);

if D.options.fitaverage
    [m_model, se_model, channels] = averageMultipleFits(cScale_model, D.channels, @mean);
else
    [m_model, se_model] = averageWithinArea(cScale_model, group_prob, @mean, 1000);
end

%% Plot
subplot('position', posc); cla; hold on
x = 1:height(channels);
tde_plotPoints(m_model(1,:)'./m_model(1,1), squeeze(se_model(1,:,:))./m_model(1,1), x, 'ci', 0, [], 40, 'r');
tde_plotPoints(m_data(1,:)'./m_data(1,1), squeeze(se_data(1,:,:))./m_data(1,1), x, 'errbar', 0, [], 40, 'k');
%tde_plotPoints(m_model(1,:)', squeeze(se_model(1,:,:)), x, 'ci', 0, [], 40, 'r');
%tde_plotPoints(m_data(1,:)', squeeze(se_data(1,:,:)), x, 'errbar', 0, [], 40, 'k');

xlim([0 max(x)+1]); ylim([0 1.5]);

set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
ylabel('residuals regression (summed)');
xlabel('Visual area');
legend('model', 'data', 'Location','NorthEast');legend boxoff

%%
subplot('position', posd); cla; hold on
x = 1:height(channels);
shift = m_model(2,1)-m_data(2,1);
%tde_plotPoints(m_model(2,:)'./max(m_model(2,:)), squeeze(se_model(2,:,:))./max(m_data(2,:)), x, 'ci', 0, [], 40, 'r');
%tde_plotPoints(m_data(2,:)'./max(m_data(2,:)), squeeze(se_data(2,:,:))./max(m_data(2,:)), x, 'errbar', 0, [], 40, 'k');

tde_plotPoints(m_model(2,:)'-shift, squeeze(se_model(2,:,:))-shift, x, 'ci', 0, [], 40, 'r');
tde_plotPoints(m_data(2,:)', squeeze(se_data(2,:,:)), x, 'errbar', 0, [], 40, 'k');
xlim([0 max(x)+1]); ylim([-0.1 1]);

set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
ylabel('R2 regression');
xlabel('Visual area');


set(findall(gcf,'-property','FontSize'),'FontSize',20)


