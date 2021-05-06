% tde_mkFigure5

% contrast response function + second pulse adaptation 
% next to
% DN model prediction + DN cascade model prediction

modelfun = @LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datatype = 'electrodeaverages_20210208T220233';
[d1] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% modelfun = @DN;
% xvalmode = 0;
% datatype = 'electrodeaverages';
% [d1] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

% modelfun = @DNCASCADE2;
% xvalmode = 0;
% datatype = 'electrodeaverages';
% [d2] = tde_loadDataForFigure(modelfun, xvalmode, datatype);
% 
% modelfun = @DNCASCADE3;
% xvalmode = 0;
% datatype = 'electrodeaverages';
% [d3] = tde_loadDataForFigure(modelfun, xvalmode, datatype);
% 
% modelfun = @DNCASCADE4;
% xvalmode = 0;
% datatype = 'electrodeaverages';
% [d4] = tde_loadDataForFigure(modelfun, xvalmode, datatype);
% 
% modelfun = @DNCASCADE5;
% xvalmode = 0;
% datatype = 'electrodeaverages';
% [d5] = tde_loadDataForFigure(modelfun, xvalmode, datatype);

srate = d1.channels.sampling_frequency(1);

%figure(1); clf
%set(gcf, 'position',  get(0, 'screensize'));
%chan_idx = 1;
saveDir = fullfile(analysisRootPath, 'figures', 'CRFISIcascade');
saveFigure = 1;

for chan_idx = 1:9
    
    figure(chan_idx); clf
    set(gcf, 'position',  get(0, 'screensize'));

    %% contrast
    conditionsOfInterest = {'CRF'};
    timepointsOfInterest = [0 0.5];

    stim_idx = find(contains(d1.stim_info.name, conditionsOfInterest));
    nStim = length(stim_idx);
    t_idx    = d1.t>timepointsOfInterest(1) & d1.t<=timepointsOfInterest(2);
    
    colors = flipud(gray(nStim+1));
    subplot(2,2,1); hold on
    for ii = 1:nStim
        plot(d1.data(t_idx,stim_idx(ii),chan_idx),'Color',colors(ii+1,:),'LineWidth', 2);
    end
    axis tight
    if chan_idx < 3 , set(gca,'ylim', [-1 20]), elseif chan_idx < 5 , set(gca,'ylim', [-1 15]), else, set(gca,'ylim', [-1 5]); end

    set(gca, 'xlim', [0 150]);
    title(sprintf('%s %s %s', d1.channels.name{chan_idx}, 'data', 'CRF'), 'fontsize', 20);
    xlabel('Time (ms)');
    
    %colors(:,[2 3]) = 0;
	subplot(2,2,2); hold on
    for ii = 1:nStim
        plot(d1.pred(t_idx,stim_idx(ii),chan_idx),'Color',colors(ii+1,:), 'LineWidth', 2);
    end
    axis tight
	title(sprintf('%s %s %s', d1.channels.name{chan_idx}, 'model', 'CRF'), 'fontsize', 20);
    xlabel('Time (ms)');
    if chan_idx < 3 , set(gca,'ylim', [-1 20]), elseif chan_idx < 5 , set(gca,'ylim', [-1 15]), else, set(gca,'ylim', [-1 5]); end

    set(gca, 'xlim', [0 150]);
    
    
    %% recovery  
    [~, ts1,w] = tde_computeISIrecovery(d1.data(:,:,chan_idx),d1.t,d1.stim_info,srate,0.5);
    [~, tspred1] = tde_computeISIrecovery(d1.pred(:,:,chan_idx),d1.t,d1.stim_info,srate,0.5);
   
    conditionsOfInterest = {'ONEPULSE4','ONEPULSE-5', 'TWOPULSE'};
    stim_idx = find(contains(d1.stim_info.name, conditionsOfInterest));
    nStim = length(find(stim_idx));

    colors = flipud(gray(nStim+1));
    subplot(2,2,3); hold on
    plot(ts1(:,1), 'k','LineWidth', 2);
    for ii = 2:nStim
        plot(ts1(:,ii),'Color', colors(ii+1,:),'LineWidth', 2);
    end
    axis tight
    if chan_idx < 3 , set(gca,'ylim', [-1 20]), elseif chan_idx < 5 , set(gca,'ylim', [-1 15]), else, set(gca,'ylim', [-1 5]); end
    set(gca, 'xlim', [0 150]);
	title(sprintf('%s %s %s', channels.name{chan_idx}, 'data', 'ISI'), 'fontsize', 20);
    xlabel('Time (ms)');
    
    %colors(:,[2 3]) = 0;
    subplot(2,2,4); hold on
    plot(tspred1(:,1), 'k','LineWidth', 2);
     for ii = 2:nStim
        plot(tspred1(:,ii),'Color', colors(ii+1,:),'LineWidth', 2);
     end
     axis tight
  	if chan_idx < 3 , set(gca,'ylim', [-1 20]), elseif chan_idx < 5 , set(gca,'ylim', [-1 15]), else, set(gca,'ylim', [-1 5]); end
    set(gca, 'xlim', [0 150]);
	title(sprintf('%s %s %s', channels.name{chan_idx}, 'model', 'ISI'), 'fontsize', 20);
    xlabel('Time (ms)');
    
    if saveFigure
        saveName = sprintf('X_%d_%s', chan_idx, channels.name{chan_idx});
        saveas(chan_idx, fullfile(saveDir, saveName), 'png');
    end
end
%
colors = flipud(turbo(5));

for chan_idx = 1:9
    
    figure(chan_idx); clf
    set(gcf, 'position',  get(0, 'screensize'));

    %% contrast
    conditionsOfInterest = {'CRF'};
    timepointsOfInterest = [-0.1 1];

    stim_idx = contains(d1.stim_info.name, conditionsOfInterest);
    t_idx    = d1.t>timepointsOfInterest(1) & d1.t<=timepointsOfInterest(2);

    subplot(2,1,1); hold on
    plot(flatten(d2.data(t_idx,stim_idx,chan_idx)),'k','LineWidth', 4);
    plot(flatten(d1.pred(t_idx,stim_idx,chan_idx)),'r', 'LineWidth', 2);
    %plot(flatten(d2.pred(t_idx,stim_idx,chan_idx)),'Color',colors(1,:),'LineWidth', 2);
    %plot(flatten(d3.pred(t_idx,stim_idx,chan_idx)),'Color',colors(2,:),'LineWidth', 2);
    %plot(flatten(d4.pred(t_idx,stim_idx,chan_idx)),'Color',colors(3,:),'LineWidth', 2);
    plot(flatten(d5.pred(t_idx,stim_idx,chan_idx)),'Color',colors(4,:),'LineWidth', 2);
    %plot(flatten(d6.pred(t_idx,stim_idx,chan_idx)),'Color',colors(5,:),'LineWidth', 2);

    nT     = length(find(t_idx));
    nStim  = length(find(stim_idx));
    axis tight
    set(gca, 'xtick',1:nT:nStim*nT, 'xticklabel', [d1.stim_info.contrast(stim_idx)*100], 'xgrid', 'on'); 
    xlabel('Contrast (%)');
    %legend({'data', 'DN', 'DNCASCADE2', 'DNCASCADE3', 'DNCASCADE4', 'DNCASCADE5'}, 'fontsize', 18, 'location', 'northwest');
    legend({'data', 'DN', 'DNCASCADE5'}, 'fontsize', 18, 'location', 'northwest');
    legend('boxoff');
    title(channels.name(chan_idx), 'fontsize', 20);
    
    %% recovery
    conditionsOfInterest = {'ONEPULSE-5', 'TWOPULSE'};
    stim_idx = contains(d1.stim_info.name, conditionsOfInterest);

    [~, ts1,w] = tde_computeISIrecovery(d2.data(:,:,chan_idx),d1.t,d1.stim_info,srate,0.5);

    [~, tspred1] = tde_computeISIrecovery(d1.pred(:,:,chan_idx),d1.t,d1.stim_info,srate,0.5);
    [~, tspred2] = tde_computeISIrecovery(d2.pred(:,:,chan_idx),d2.t,d2.stim_info,srate,0.5);
    [~, tspred3] = tde_computeISIrecovery(d3.pred(:,:,chan_idx),d3.t,d3.stim_info,srate,0.5);
    [~, tspred4] = tde_computeISIrecovery(d4.pred(:,:,chan_idx),d4.t,d4.stim_info,srate,0.5);
    [~, tspred5] = tde_computeISIrecovery(d5.pred(:,:,chan_idx),d5.t,d5.stim_info,srate,0.5);
    %[~, tspred6] = tde_computeISIrecovery(d6.pred(:,:,chan_idx),d6.t,d6.stim_info,srate,0.5);

    subplot(2,1,2); hold on

    plot(flatten(ts1),'k','LineWidth', 4);
    plot(flatten(tspred1),'r', 'LineWidth', 2);
    %plot(flatten(tspred2),'Color',colors(1,:),'LineWidth', 2);
    %plot(flatten(tspred3),'Color',colors(2,:),'LineWidth', 2);
    %plot(flatten(tspred4),'Color',colors(3,:),'LineWidth', 2);
    plot(flatten(tspred5),'Color',colors(4,:),'LineWidth', 2);
    %plot(flatten(tspred6),'Color',colors(5,:),'LineWidth', 2);

    axis tight

    nT     = w/(1/srate);
    nStim  = length(find(stim_idx))+1;
    set(gca, 'xtick',1:nT:nStim*nT, 'xticklabel', [0; d1.stim_info.ISI(stim_idx)*1000], 'xgrid', 'on'); 

    xlabel('ISI (ms)');
    
    if saveFigure
        saveName = sprintf('A_%d_%s', chan_idx, channels.name{chan_idx});
        saveas(chan_idx, fullfile(saveDir, saveName), 'png');
    end
end


%%
modelfun = @LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datatype = 'electrodeaverages';
[d1] = tde_loadDataForFigure(modelfun, xvalmode, datatype);


