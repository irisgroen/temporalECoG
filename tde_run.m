% script 

tic

% load or (re)compute the processed data
reComputeFlag = false; 
[data] = tde_getData(reComputeFlag);

% select epochs and channels, average trials within stimulus condition 
opts = [];
opts.doplots         = false;
opts.average_trials  = true;
opts.normalize_data  = false;
opts.average_elecs   = true;
opts.sort_channels   = true;
[data2fit, channels, stimnames, t, srate] = tde_selectData(data, [], opts);

% generate stimulus timecourses
[stim_ts] = tde_generateStimulusTimecourses(stimnames,t);

toc

% TEMPORARY: plot final selected data
figure, sz = ceil(sqrt(size(data2fit,3)));
for ii = 1:size(data2fit,3)
    subplot(sz,sz,ii); plot(t,data2fit(:,:,ii), 'LineWidth', 2);
    if opts.average_elecs
        title(sprintf('%s (n = %d)', ...
            channels.name{ii}, channels.number_of_elecs{ii})); 
    else
        title(sprintf('%s %s %s %s', ...
            channels.bensonarea{ii}, channels.wangarea{ii}, channels.subject_name{ii}, channels.name{ii})); 
    end
    yaxlims = get(gca, 'YLim');
    line([0 0], [0 yaxlims(2)], 'Color', 'k', 'LineStyle', ':')
    %line([t(1) t(end)], [0 0],'Color', 'k', 'LineStyle', ':');
    set(gca, 'YLim', [0 yaxlims(2)]);
end
set(gcf, 'Position', [400 200 2000 1200]);

%% fitting

% define electrode (temporary)
%ele = 54; % 58; % 2;%54; 
%smallData = data2fit(:,:,ele);

% define model
modelfuns = tde_modelTypes();
modelfun = modelfuns{4}; % eg 'DNCASCADE'; 

tic
[params, pred] = tde_fitModel(modelfun, data2fit, stim_ts, srate);
toc

%% evaluation of fits

[results] = tde_evaluateModelFit(modelfun, params, data2fit, pred);

%% TEMPRARY plotting

% % visualization 1
% figure;
% subplot(2,1,1);plot(t,smallData, 'LineWidth', 3); title('data')
% subplot(2,1,2);plot(t,pred, 'LineWidth', 3); title('model fit')

% visualization 2
conditionsOfInterest = {'CRF', 'ONEPULSE', 'TWOPULSE'};
nCond = length(conditionsOfInterest);
for kk = 1:size(data2fit,3)
    figure('Name', channels.name{kk});
    d = data2fit(:,:,kk);
    maxresp = max(d(:));
    p = pred(:,:,kk);
    for ii = 1:length(conditionsOfInterest)
        inx = contains(stimnames, conditionsOfInterest{ii});
        subplot(nCond,1,ii); hold on
        h = plot(flatten(stim_ts(:,inx))*maxresp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        plot(flatten(d(:,inx)), 'k', 'LineWidth', 3); 
        plot(flatten(p(:,inx)), 'r', 'LineWidth', 3);
        legend('data', 'model fit'); 
        title(sprintf('%s r2 = %0.2f',conditionsOfInterest{ii}, mean(results.rSquare(inx,kk))));
        axis tight   
        set(gca, 'XTick',1:length(t):length(find(inx))*length(t), 'XTickLabel', []);
    end
    set(gcf, 'Position', [400 200 1000 1200]);
end

%% plot summaries

figure;hold on
derivedTitles = {'time2peak', 'R_asymp'};
fittedTitles = {'tau1', 'weight','tau2', 'n', 'sigma'};
nChans = height(channels);
subplot(2,5,1); plot(1:nChans,mean(results.rSquare), '.b', 'MarkerSize', 50, 'LineStyle', 'none')
set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
title('explained variance'); xlabel('visual area');  ylabel('R2'); set(gca, 'fontsize', 16);

for p = 1:2
    subplot(2,5,p+1);
    plot(1:nChans,results.derivedPrm(p,:), '.r', 'MarkerSize', 50, 'LineStyle', 'none')
    set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
    title(derivedTitles{p}); xlabel('visual area');  ylabel('parameter value'); set(gca, 'fontsize', 16);
end
for p = 1:5
    subplot(2,5,p+5);
    plot(1:size(data2fit,3),results.fittedPrm(p,:), '.k', 'MarkerSize', 50, 'LineStyle', 'none')
    set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
    title(fittedTitles{p}); xlabel('visual area'); ylabel('parameter value');set(gca, 'fontsize', 16);
end
set(gcf, 'Position', [400 200 2000 1200]);

%%
% -- which models?
% ----- DN (flavors: uniphasic, biphasic, fixed exponent or not)
% ----- DN cascade?
% ----- Two temporal channels (flavors: HH, Stigliani 1, Stigliani 2)
% ----- DN-like models:
% --------- Heeger 1993

% next step, tde_analyzeModelFits
% input: model fits 
% ---> put deriveParams and R2 in here (rather than in tde_modelFit)
% outputs: plots, stats

% function tde_plotModelFits 


% How about the PRFs --> separate pipeline (like this one), read in fits
% from file (e.g. add to channel table), prf_getData, prf_selectData,
% prf_fitModel



