%% 1: Load the ECoG data and stimulus description

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

tde_plotData(data2fit, channels, t, opts)

%% 2. Model fitting

% define electrode (temporary)
ele = 59; % 58; % 2;%54; 
%smallData = data2fit(:,:,ele);
smallData = data2fit;
% define model
modelfuns = tde_modelTypes();

tic
modelfun = modelfuns{5}; % TTCSTIG19
[params1, pred1] = tde_fitModel(modelfun, smallData, stim_ts, srate);
toc
tic
modelfun = modelfuns{1}; % DN
[params2, pred2] = tde_fitModel(modelfun, smallData, stim_ts, srate);
toc
%% evaluation of fits

[results1] = tde_evaluateModelFit(modelfun, params1, smallData, pred1);
[results2] = tde_evaluateModelFit(modelfun, params2, smallData, pred2);
[results] = tde_evaluateModelFit(modelfun, params2, smallData, pred2);

%% TEMPORARY plotting

% % visualization 1
% figure;
% subplot(2,1,1);plot(t,smallData, 'LineWidth', 3); title('data')
% subplot(2,1,2);plot(t,pred, 'LineWidth', 3); title('model fit')

% visualization 2
conditionsOfInterest = {'CRF', 'ONEPULSE', 'TWOPULSE'};
nCond = length(conditionsOfInterest);
for kk = 1:size(smallData,3)
    figure('Name', channels.name{kk});
    d = smallData(:,:,kk);
    maxresp = max(d(:));
    p = pred1(:,:,kk);
    p1 = pred2(:,:,kk);
    for ii = 1:length(conditionsOfInterest)
        inx = contains(stimnames, conditionsOfInterest{ii});
        subplot(nCond,1,ii); hold on
        h = plot(flatten(stim_ts(:,inx))*maxresp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        plot(flatten(d(:,inx)), 'k', 'LineWidth', 2); 
        plot(flatten(p(:,inx)), 'b', 'LineWidth', 2);
        plot(flatten(p1(:,inx)), 'r', 'LineWidth', 2);
        legend('data', 'TTCSTIG19','DN'); 
        title(sprintf('%s r2 DN = %0.2f , r2 TTCSTIG19 = %0.2f',conditionsOfInterest{ii}, mean(results2.rSquare(inx,kk)), mean(results1.rSquare(inx,kk))));
        axis tight   
        set(gca, 'XTick',1:length(t):length(find(inx))*length(t), 'XTickLabel', []);
    end
    set(gcf, 'Position', [400 200 1800 1200]);
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

% function tde_plotModelFits 

% How about the PRFs --> separate pipeline (like this one), read in fits
% from file (e.g. add to channel table), prf_getData, prf_selectData,
% prf_fitModel



