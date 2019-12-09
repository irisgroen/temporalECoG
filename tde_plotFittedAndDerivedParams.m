function tde_plotFittedAndDerivedParams(results, channels)

nChans = height(channels);
nModels     = size(results,2);
nDatasets   = size(results(1).params,2);

%% Plot prediction for sustained stimulus
for ii = 1:nDatasets
    figure('Name', sprintf('%s %s', channels.name{ii}, 'derived predictions')); hold on
    l = cell(1,nModels);
    for kk = 1:nModels
        plot(results(kk).derivedPred(:,ii), 'LineWidth', 2);
        set(gca, 'Xlim', [0 1000]);
        l{kk} = sprintf('%s t2p = %0.2f rasymp = %0.2f', func2str(results(kk).model), results(kk).derivedPrm(1,ii), results(kk).derivedPrm(2,ii));
    end
	legend(l); 
    set(gca, 'FontSize', 16);
	title('Model predictions for sustained stimulus');
end

%% plot fitted and derived Params

% use this to determine if data was averaged across elecs prior to fit
if ~isfield(summary(channels), 'number_of_elecs')
    %plotErrorBar = 0;
    % to do: average params across electrodes from same visual area
end


for kk = 1:nModels
    figure('Name', sprintf('%s %s', func2str(results(kk).model), 'parameters')); hold on
    derivedTitles = {'time2peak', 'Rasymp'};
    %fittedTitles = {'tau1', 'weight','tau2', 'n', 'sigma'};
    tmp = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', func2str(results(kk).model))));
    paramNames = strsplit(tmp.params,',');
    nParams = length(paramNames);
    
    subplot(2,nParams,1); plot(1:nChans,results(kk).rSquareConc, '.b', 'MarkerSize', 50, 'LineStyle', 'none')
    set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
    title('explained variance'); xlabel('visual area');  ylabel('R2'); set(gca, 'fontsize', 16);

    for p = 1:2
        subplot(2,nParams,p+1);
        plot(1:nChans,results(kk).derivedPrm(p,:), '.r', 'MarkerSize', 50, 'LineStyle', 'none')
        set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
        title(derivedTitles{p}); xlabel('visual area');  ylabel('parameter value'); set(gca, 'fontsize', 16);
    end
    for p = 1:nParams
        subplot(2,nParams,p+nParams);
        plot(1:nChans,results(kk).params(p,:), '.k', 'MarkerSize', 50, 'LineStyle', 'none')
        set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
        title(paramNames{p}); xlabel('visual area'); ylabel('parameter value');set(gca, 'fontsize', 16);
    end
    set(gcf, 'Position', [400 200 2000 1200]);

end


end