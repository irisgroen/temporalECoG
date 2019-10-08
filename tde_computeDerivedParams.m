function derivedPrm = tde_computeDerivedParams(objFunction, prm)
% Simulate model response to a long sustained stimulus
% Extract two summary statistics:
% - t2pk = time to peak 
% - r_asymp = magnitude or response at asymptote

%% compute derived parameters

t    = 0.001 : 0.001 : 10;
stim = ones(length(t),1);
srate = 1/median(diff(t)); % samples per second

%prm(end) = 1; % Jing's code appears to fix the gain to 1, do we want this?
[~, pred] = objFunction(prm, [], stim, srate);

derivedPrm.t2pk = t(pred == max(pred));
derivedPrm.r_asymp = pred(end);

% % compute model responses
% for k = 1 : size(prm, 1)
%     
%    [~, rsp(k, :)] = dn_computeFineFit(prm(k, :), stim, stim, t, irfType);
%    
%    % compute time to peak
%    derivedPrm.t2pk(k) = find(rsp(k, :) == max(rsp(k, :)));
%    
%    % compute asymptotic response
%    derivedPrm.r_asymp(k) = rsp(k, end);
% end

end