paramNames = {'tau1','tau2','n_irf','weight','shift','scale','n','sigma','tau_a'};


figure;hold on
for ii = 1:size(params,1)
    subplot(3,3,ii);
    plot(params(ii,:), 'LineWidth', 2);
    title(paramNames{ii});
end
    