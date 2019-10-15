% script 

% note: sub-beilen not included at the moment because of inconsistencies in
% the formatting of the events.tsv files, need to sync with Gio.
tic

% load (0) or (re)compute (1)
reComputeFlag = true; 
[data] = tde_getData(reComputeFlag);

% select epochs and channels, average trials within stimulus condition 
opts = [];
opts.doplots         = true;
opts.normalize_data  = true;
opts.average_elecs   = false;
[data2fit, channels, stimnames, t] = tde_selectData(data, [], opts);

% generate stimulus timecourses
[stim_ts] = tde_generateStimulusTimecourses(stimnames,t);

srate = 1/median(diff(t)); % samples per second
toc

%% fitting

ele = 54; 
smallData = data2fit(:,:,ele);
figure;

subplot(2,2,1);plot(t,smallData); title('data')


opts = [];
opts.srate = srate;
opts.x0   = [0.03, 0, 0.07, 1.5, 0.15, 0.06, 1];
opts.lb   = [0, 0, 0, 0, 0, 0, 0];
opts.ub   = [1, 0, 1, 10, 1, 1, 1];

tic
[results1, pred1] = tde_fitModel(@DNmodel, smallData, stim_ts, opts);
toc
subplot(2,2,3);plot(t,pred1); title('tdefitmodel')


% compare with Jings original code:
tic
opts.x0   = [0.03, 0.07, 1.5, 0.15, 0.06, 1];
opts.lb   = [0, 0, 0, 0, 0, 0];
opts.ub   = [1, 1, 10, 1, 1, 1];
fprintf('[%s] Fitting dn_DNmodel \n', mfilename);
prm = fminsearchbnd(@(x) dn2_fineFitCtrstDur(x, smallData', t, stim_ts'), opts.x0, opts.lb, opts.ub);

prm_tofit = [prm(1), 0, prm(2 : end)];
pred2 = dn_DNmodel(prm_tofit, stim_ts', t);
derived_prm = dn_computeDerivedParams(prm, 'uniphasic');

results2.derivedPrm(1,1) = derived_prm.t2pk;
results2.derivedPrm(2,1) = derived_prm.r_asymp;
results2.fittedPrm = prm';
pred2 = pred2./max(pred2(:));
for k = 1:17
    results2.rSquare(k,1) = corr(pred2(k,:)',smallData(:,k)).^2;  
end
disp(mean(results2.rSquare));
toc
subplot(2,2,4);plot(t,pred2');title('dn_DNmodel')

% input: model: function handle, data, eventCodes, stimulusTimeCourses, opts
% output: model fits
% -- which models?
% ----- DN (flavors: uniphasic, biphasic, fixed exponent or not)
% ----- DN cascade?
% ----- Two temporal channels (flavors: HH, Stigliani 1, Stigliani 2)
% ----- DN-like models:
% --------- Heeger 1993

% next step, tde_analyzeModelFits
% input: model fits (can be multiple?)
% outputs: plots, stats

% function tde_plotModelFits 


% How about the PRFs --> separate pipeline (like this one), read in fits
% from file (e.g. add to channel table), prf_getData, prf_selectData,
% prf_fitModel



