%% compare with Jings original code:
tic
opts.x0   = [0.03, 0.07, 1.5, 0.15, 0.06, 1];
opts.lb   = [0, 0, 0, 0, 0, 0];
opts.ub   = [1, 1, 10, 1, 1, inf];
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
