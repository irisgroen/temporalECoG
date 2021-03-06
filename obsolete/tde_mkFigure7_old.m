% tde_mkFigure 7

modelfun = @LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;

datatype = 'electrodeaverages';
d1 = tde_loadDataForFigure(modelfun, xvalmode, datatype);

datatype = 'individualelecs';
d2 = tde_loadDataForFigure(modelfun, xvalmode, datatype);

[results] = tde_evaluateModelFit(d2);

[stim_ts, stim_info] = tde_generateStimulusTimecourses(d1.options.stimnames,d1.t);

%% Plot specs
% Subplot positions: % [left bottom width height]

posa1 = [0.05  0.75 0.9 0.22];
posa2 = [0.05  0.50 0.9 0.22];
posb =  [0.05  0.1 0.18 0.3];
posc =  [0.28  0.1 0.18 0.3];
posd =  [0.53  0.1 0.18 0.3];
pose =  [0.77  0.1 0.18 0.3];

% posb =  [0.05  0.1 0.25 0.3];
% posc =  [0.375 0.1 0.25 0.3];
% posd =  [0.7   0.1 0.25 0.3];

figure(1); clf
set(gcf, 'position',  get(0, 'screensize'));

% Panel A: time courses comparison V1 and V3b

conditionsOfInterest = {'CRF'};
timepointsOfInterest = [-0.1 1.2];
areasOfInterest      = {'V1', 'V3b'};

stim_idx  = contains(stim_info.name, conditionsOfInterest);
t_idx     = d1.t>timepointsOfInterest(1) & d1.t<=timepointsOfInterest(2);
chan_idx1 = find(strcmp(d1.channels.name, areasOfInterest{1}));
chan_idx2 = find(strcmp(d1.channels.name, areasOfInterest{2}));

% Top: Data

subplot('position', posa1); hold on

s = stim_ts(t_idx,stim_idx);
c1 = d1.data(t_idx,stim_idx,chan_idx1);
c2 = d1.data(t_idx,stim_idx,chan_idx2);
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

ylabel('Response (normalized)'); %title('Broadband responses to increasing durations'); 
legend(areasOfInterest, 'location', 'northeast');
legend('boxoff');

% Bottom: Model

subplot('position', posa2); hold on

s = stim_ts(t_idx,stim_idx);
c1 = d1.pred(t_idx,stim_idx,chan_idx1);
c2 = d1.pred(t_idx,stim_idx,chan_idx2);
c1 = c1/max(c1(:,end));
c2 = c2/max(c2(:,end));

hs = plot(s(:), 'color', [0.7 0.7 0.7], 'linewidth', 1);
hd = plot(c1(:), 'r-', 'linewidth', 2);
hp = plot(c2(:), 'Color', [1 0.5 0.5], 'linewidth', 2);
set(get(get(hs,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

set(gca, 'xtick',1:size(c1,1):size(c1,2)*size(c1,1), 'xticklabel', stim_info.ISI(stim_idx)*1000);
set(gca, 'ytick', [0 1]);
box off,  axis tight

ylim([-0.2 1.5]);
xlim([-20 length(s(:)) + 20]);

xlabel('Stimulus duration (ms)'); %ylabel('Response (normalized)'); %title('Broadband responses to increasing durations'); 
legend(areasOfInterest, 'location', 'northeast');
legend('boxoff');

%% Panel B:recovery for different areas superimposed
[~, tspred1,w] = tde_computeISIrecovery(d1.data,d1.t,d1.stim_info,d1.srate,0.5, [], 'max');
[~, tspred2] = tde_computeISIrecovery(d1.pred,d1.t,d1.stim_info,d1.srate,0.5, [], 'max');
COI = 8;

cmap1      = brewermap(height(d1.channels)+2, '*Greys');
cmap1      = cmap1(1:height(d1.channels),:);
cmap2      = brewermap(height(d1.channels)+2, '*OrRd');
cmap2      = cmap2(1:height(d1.channels),:);

t0 = d1.t(find(d1.t>0,1));
x = t0:(1/d1.srate):w;
x = floor(x*1000)/1000; 
s_off = find(x>0.133,1);
s = zeros(size(tspred1,1),1);
s(2:s_off) = 1;

subplot('position', posb); hold on
hs = plot(x,s, 'color', [0.7 0.7 0.7], 'linewidth', 1);
hs.Annotation.LegendInformation.IconDisplayStyle = 'off';

ts2 = squeeze(tspred2(:,COI,:));
p2 = plot(x, ts2, 'LineWidth',2);
set(p2, {'color'}, num2cell(cmap2,2));

xlabel('Time (s)'); ylabel('Response (normalized)');
legend(d1.channels.name);
xlim([-0.1 w+0.10]); 
%ylim([0 1.1]);
%set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
%set(gca, 'ytick', [0 0.1 0.2 0.3]);
legend('boxoff')

subplot('position', posc); hold on
hs = plot(x,s, 'color', [0.7 0.7 0.7], 'linewidth', 1);
hs.Annotation.LegendInformation.IconDisplayStyle = 'off';

ts2 = squeeze(tspred2(:,COI,:)./max(tspred1(:,1,:)));
p2 = plot(x, ts2, 'LineWidth',2);
set(p2, {'color'}, num2cell(cmap2,2));

xlabel('Time (s)'); ylabel('Response (normalized)');
xlim([-0.1 w+0.10]); 

subplot('position', posd); hold on
hs = plot(x,s, 'color', [0.7 0.7 0.7], 'linewidth', 1);
hs.Annotation.LegendInformation.IconDisplayStyle = 'off';

ts2 = squeeze(tspred2(:,COI,:));
ts2 = ts2./max(ts2);
p2 = plot(x, ts2, 'LineWidth',2);
set(p2, {'color'}, num2cell(cmap2,2));

xlabel('Time (s)'); ylabel('Response (normalized)');
xlim([-0.1 w+0.10]); 

set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% Panel D, C, E: derived parameters summarized across areas
[~, channels, group_prob] = groupElecsByVisualArea(d2.channels, 'probabilisticresample');   
[m, se] = averageWithinArea(results.derived.params, group_prob, [], 10000);

subplot('position', posb); hold on
x = 1:height(channels);
tde_plotPoints(m(4,:)', squeeze(se(4,:,:)), x, 'errbar', 0)
xlim([0 max(x)+1]); ylim([0 0.8]);
set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
set(gca, 'ytick', 0:0.2:0.8);
ylabel('Time to recover 50% (s)');
xlabel('Visual area');

subplot('position', posc); hold on
x = 1:height(channels);
tde_plotPoints(m(5,:)', squeeze(se(4,:,:)), x, 'errbar', 0)
xlim([0 max(x)+1]); ylim([0 0.8]);
set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
set(gca, 'ytick', 0:0.2:0.8);
ylabel('Time to recover 80% (s)');
xlabel('Visual area');

subplot('position', posd); hold on
x = 1:height(channels);
tde_plotPoints(m(6,:)', squeeze(se(4,:,:)), x, 'errbar', 0)
xlim([0 max(x)+1]); ylim([0 0.8]);
set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
set(gca, 'ytick', 0:0.2:0.8);
ylabel('Time to recover 100 %(s)');
xlabel('Visual area');

set(findall(gcf,'-property','FontSize'),'FontSize',20)

%% Extended data: individual subjects
figure(2); clf
set(gcf, 'position',  get(0, 'screensize'));

x = 1:height(channels);

ylabels = {'Time to peak', 'Ratio sustained/transient', 'Fwhm', 'Time to recover'};

xlabel('Visual area');

subjectNames = unique(d2.channels.subject_name);

cmap = brewermap(length(subjectNames),'Set2');
%cmap(:,[2 3]) = 0.5;

nRows = length(ylabels);
nCols = length(subjectNames);
minY = min(results.derived.params, [], 2);
maxY = max(results.derived.params, [], 2);

c = 0;
for jj = 1:nRows
    
    for ii = 1:nCols
        
        c = c+1;
        idx = contains(d2.channels.subject_name, subjectNames{ii});
        [~, ~, group_prob] = groupElecsByVisualArea(d2.channels(idx,:), 'probabilisticresample');  
        [m, se] = averageWithinArea(results.derived.params(:,idx), group_prob, [], 1000);
        
        subplot(nRows,nCols,c);hold on

        ix = find(~isnan(m(1,:)));
        [hp, hc] = tde_plotPoints(m(jj,ix)', squeeze(se(jj,ix,:)), x(ix), 'errbar', 0, [], 20);
        hc.FaceColor = cmap(ii,:);
        hp.Color = cmap(ii,:);
        %hp.Color = cmap(ii,:);
        %hp.LineWidth = 1;
        xlim([0 max(x)+1]); 
        %ylim([minY(jj)-0.1*minY(jj) maxY(jj)+0.1*maxY(jj)]);
        if jj == 1, ylim([0 0.4]), end
        if jj == 2, ylim([0 0.7]), end
        if jj == 3, ylim([0 0.3]), end
        if jj == 4, ylim([0 1.5]), end

        %set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
        if ii > 1, set(gca, 'ytick', []); end
        if ii == 1, ylabel(ylabels{jj}); end

        if jj == nRows
            set(gca, 'xtick', x, 'xticklabel', channels.name, 'xticklabelrotation', 45);
        else
            set(gca, 'xtick', []); 
        end
        if jj == 1
            title(subjectNames{ii});
        end

    end

%     [~, channels, group_prob] = groupElecsByVisualArea(d2.channels, 'probabilisticresample');   
%     [m, se] = averageWithinArea(results.derived.params, group_prob, [], 1000);
%     h = tde_plotPoints(m(jj,:)', squeeze(se(jj,:,:)), x, 'errbar', 0, '-');
%     h.LineWidth = 2;
%     h.MarkerSize = 50;
%     ylabel(ylabels{jj});
%     if jj == 1
%         legend([subjectNames; 'all']);
%     end


  
end
set(findall(gcf,'-property','FontSize'),'FontSize',14)
