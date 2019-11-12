% script 

tic

% load or (re)compute the processed data
reComputeFlag = false; 
[data] = tde_getData(reComputeFlag);

% select epochs and channels, average trials within stimulus condition 
data_opts = [];
data_opts.doplots         = false;
data_opts.normalize_data  = false;
data_opts.average_elecs   = false;
[data2fit, channels, stimnames, t, srate] = tde_selectData(data, [], data_opts);

% generate stimulus timecourses
[stim_ts] = tde_generateStimulusTimecourses(stimnames,t);

toc

% debug: check data time series
% figure, sz = ceil(sqrt(size(data2fit,3)));
% for ii = 1:size(data2fit,3)
%     subplot(sz,sz,ii); plot(data2fit(:,17,ii), 'LineWidth', 3); 
%     title(ii); 
% end

%% fitting

% define electrode (temporary)
ele = 54; % 58; % 2;%54; 
smallData = data2fit(:,:,ele);

% define model
modelType = 'DNCASCADE'; 

modelfun = str2func(sprintf('%smodel', modelType));

tic
[results, pred] = tde_fitModel(modelfun, smallData, stim_ts, srate);
toc

% plotting

% visualization 1
figure;
subplot(2,1,1);plot(t,smallData, 'LineWidth', 3); title('data')
subplot(2,1,2);plot(t,pred, 'LineWidth', 3); title('model fit')

% visualization 2
conditionsOfInterest = {'CRF','ONEPULSE', 'TWOPULSE'};
nCond = length(conditionsOfInterest);
figure;
for ii = 1:length(conditionsOfInterest)
    inx = contains(stimnames, conditionsOfInterest{ii});
    subplot(nCond,1,ii); hold on
    plot(flatten(smallData(:,inx)), 'k', 'LineWidth', 3); 
    plot(flatten(pred(:,inx)), 'r', 'LineWidth', 3);
    legend('data', 'model fit'); 
    title(conditionsOfInterest{ii});
    axis tight   
    set(gca, 'XTick',1:length(t):length(find(inx))*length(t), 'XTickLabel', []);
end

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



