function diff = dn2_fineFitCtrstDur(prm, data, t, stim)
% [prm, prd, r2] = dn2_fineFitCtrstDur(prm, data, t, stim, seed)
%
% INPUTS ------------------------------------------------------------------
% prm : [tau1, weight, tau2, n, sigma, shift, scale]
% t   : a vector in unit of ms
% stim: num.stim x time course

% EXAMPLE -----------------------------------------------------------------
% prm = [0.07, 0, 0.1, 2, 0.05, 0.01, 1];

%% PRE-DEFINED VARIABLES


%% GENERATE DN MODEL PREDICTION

prm_tofit = [prm(1), 0, prm(2 : end)];

pred = dn_DNmodel(prm_tofit, stim, t);
pred = pred./max(pred(:));

diff = sum((pred(:) - data(:)).^2);

%% VISUALIZE

% figure (100), clf
% subplot(1, 2, 1), imagesc(pred), 
% subplot(1, 2, 2), imagesc(data), drawnow

% figure (100), clf
% for k = 1 : 17
%    subplot(17, 1, k)
%    plot(data(k, :)), hold on
%    plot(pred(k, :)), axis tight, box off, ylim([0, 1]), drawnow
% end

end