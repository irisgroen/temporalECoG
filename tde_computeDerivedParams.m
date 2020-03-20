function [derivedPrm, pred, t] = tde_computeDerivedParams(objFunction, prm)
% Simulate model response to a long sustained stimulus
% Extract two summary statistics:
% - t2pk = time to peak 
% - r_asymp = magnitude or response at asymptote

%% Compute derived parameters Time2Peak and Rasymptote

t    = 0.001 : 0.001 : 10;
stim = ones(length(t),1);
srate = 1/median(diff(t)); % samples per second

[~, pred] = objFunction(prm, [], stim, srate);

[~,x] = max(pred);
derivedPrm.t2pk    = t(x);
derivedPrm.r_asymp = pred(end)/max(pred);

%% Compute derived parameter Rdouble

T = 3;
stim = zeros(T*1000,2);
stim(1:20,1) = 1;
stim(1:40,2) = 1;

[~, pred] = objFunction(prm, [], stim, srate);
rsp = sum(pred, 1);
derivedPrm.r_double = rsp(2)/(2*rsp(1));

% % debug
% figure;plot(stim, 'LineWidth', 2)
% %hold on; plot(pred./max(pred), 'LineWidth', 2)
% hold on; plot(pred, 'LineWidth', 2)
% %ylim([-1 2])
%% Compute derived parameter Tisi

T = 5;
stim = zeros(T*1000,1);
stim(1 : 100) = 1;

thresh = 2-1/exp(1);

% compute linear prediction
[~, pred] = objFunction(prm, [], stim, srate);
linpred = sum(pred,1);

% compute t_isi
options = gaoptimset('Display', 'off');
    
[t_isi, ~, exitFlg] = ga(@(x) fit_t_isi(x, prm, thresh, stim, objFunction, srate, linpred), ...
    1, [], [], [], [], 1, 4500, [], 1,  options);
if exitFlg == 0, t_isi = nan; end

derivedPrm.t_isi = t_isi;

% % debug
% figure;
% subplot(2,1,1);hold on
% plot(stim, 'LineWidth', 2)
% plot(pred, 'LineWidth', 2)
% set(gca, 'Xlim', [0 1000]);
% stim_on  = find(stim == 1);
% stim_lth = length(stim_on);
% stim_end = stim_on(end);
% pulse2_st  = stim_end + t_isi + 1;
% pulse2_end = stim_end + t_isi + 1 + stim_lth;
% stim2 = stim;
% stim2(pulse2_st : pulse2_end) = 1;
% [~, pred2] = objFunction(prm, [], stim2, srate);
% subplot(2,1,2);hold on
% plot(stim2, 'LineWidth', 2)
% plot(pred2, 'LineWidth', 2)
% set(gca, 'Xlim', [0 1000]);

end

function residual = fit_t_isi(gap, prm, thresh, stim, objFunction, srate, linpred)


% create new stimulus
stim_on  = find(stim == 1);
stim_lth = length(stim_on);
stim_end = stim_on(end);

% create the second pulse
pulse2_st  = stim_end + gap + 1;
pulse2_end = stim_end + gap + 1 + stim_lth;

stim(pulse2_st : pulse2_end) = 1;

%% compute

[~, pred] = objFunction(prm, [], stim, srate);
rsp = sum(pred,1);

%% compute residual

residual = abs(linpred * thresh - rsp);

%%
% figure(1), clf
% plot(linpred*thresh, 'ro'), hold on
% plot(rsp, 'bo'), drawnow

end


