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

% plot average response per stimulus for selected data
tde_plotData(data2fit, channels, t, opts)

%% 2. Model fitting

% define electrode (temporary)
% ele = 59; % 58; % 2;%54; 
% smallData = data2fit(:,:,ele);
smallData = data2fit(:,:,1);
% define model
modelfuns = tde_modelTypes();

options = struct();
options.xvalmode = 'stimuli';
options.display  = 'final';

tic
modelfun = modelfuns{1}; % DN
[params2, pred2] = tde_fitModel(modelfun, stim_ts, smallData, srate, options);
toc
tic
modelfun = modelfuns{5}; % TTCSTIG19
[params1, pred1] = tde_fitModel(modelfun, stim_ts, smallData, srate, options);
toc


%% evaluation of fits

modelfun = modelfuns{5}; % TTCSTIG19
[results1] = tde_evaluateModelFit(modelfun, params1, smallData, pred1);
modelfun = modelfuns{1}; % DN
[results2] = tde_evaluateModelFit(modelfun, params2, smallData, pred2);

figure;hold on
plot(results1.derivedPred, 'LineWidth', 2);
plot(results2.derivedPred, 'LineWidth', 2);
set(gca, 'Xlim', [0 1000]);
leg1 = sprintf('TTCSTIG t2p = %0.2f rasymp = %0.2f', results1.derivedPrm(1), results1.derivedPrm(2));
leg2 = sprintf('DN t2p = %0.2f rasymp = %0.2f', results2.derivedPrm(1), results2.derivedPrm(2));
legend({leg1,leg2}); 
set(gca, 'FontSize', 12);
title('Model predictions for sustained stimulus');

%% TEMPORARY plotting

% visualization 2
conditionsOfInterest = {'CRF', 'ONEPULSE', 'TWOPULSE'};
nCond = length(conditionsOfInterest);
for kk = 1:size(smallData,3)
    figure('Name', channels.name{kk});
    d = smallData(:,:,kk);
    maxresp = max(d(:));
    p1 = pred1(:,:,kk);
    p2 = pred2(:,:,kk);
    for ii = 1:length(conditionsOfInterest)
        inx = contains(stimnames, conditionsOfInterest{ii});
        subplot(nCond,1,ii); hold on
        h = plot(flatten(stim_ts(:,inx))*maxresp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        plot(flatten(d(:,inx)), 'k', 'LineWidth', 2); 
        plot(flatten(p1(:,inx)), 'b', 'LineWidth', 2);
        plot(flatten(p2(:,inx)), 'r', 'LineWidth', 2);
        legend('data', 'TTCSTIG19','DN'); 
        title(sprintf('%s r2 DN = %0.2f , r2 TTCSTIG19 = %0.2f',conditionsOfInterest{ii}, mean(results2.rSquare(inx,kk)), mean(results1.rSquare(inx,kk))));
        axis tight   
        set(gca, 'XTick',1:length(t):length(find(inx))*length(t), 'XTickLabel', []);
    end
    set(gcf, 'Position', [400 200 1800 1200]);
end

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



