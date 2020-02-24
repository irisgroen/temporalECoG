function [derivedPrm, pred, t] = tde_computeDerivedParams(objFunction, prm)
% Simulate model response to a long sustained stimulus
% Extract two summary statistics:
% - t2pk = time to peak 
% - r_asymp = magnitude or response at asymptote

%% compute derived parameters

t    = 0.001 : 0.001 : 10;
stim = ones(length(t),1);
srate = 1/median(diff(t)); % samples per second

[~, pred] = objFunction(prm, [], stim, srate);

[~,x] = max(pred);
derivedPrm.t2pk    = t(x);
derivedPrm.r_asymp = pred(end)/max(pred);

end