function [results] = tde_plotFittedAndDerivedParams(results, channels)

nChans = height(channels);
nModels     = size(results,2);
nDatasets   = size(results(1).params,2);

%% Plot prediction for sustained stimulus
colors = {'r', 'b', 'c', 'm'}; % assuming we'll never plot >4 model fits at a time

figure('Name', sprintf('%s', 'Derived predictions')); 
for ii = 1:nDatasets
    subplot(ceil(sqrt(nDatasets)),ceil(sqrt(nDatasets)),ii); hold on
    l = cell(1,nModels);
    for kk = 1:nModels
        plot(results(kk).derivedPred(:,ii), 'Color', colors{kk}, 'LineWidth', 2);
        set(gca, 'Xlim', [0 1000]);
        l{kk} = sprintf('%s t2p = %0.2f rasymp = %0.2f', func2str(results(kk).model), results(kk).derivedPrm(1,ii), results(kk).derivedPrm(2,ii));
    end
	legend(l); 
    title(channels.name{ii});
    set(gca, 'FontSize', 14);
    set(gcf, 'Position', [400 200 2000 1200]);
end

%% plot parameters

% Determine if data was averaged across elecs prior to fit; if not, average
% derived and fitted parameters now
if ~isfield(summary(channels), 'number_of_elecs')
    newresults = results;
    [INX, channels] = groupElecsByVisualArea(channels);    
    nAreas = size(INX,2);
    for kk = 1:nModels
        for jj = 1:nAreas
            newresults(kk).params(:,jj) = mean(results(kk).params(:,INX{jj},2));
            newresults(kk).derivedPrm(:,jj) = mean(results(kk).derivedPrm(:,INX{jj},2));
            newresults(kk).rSquareConc(:,jj) = mean(results(kk).rSquareConc(:,INX{jj},2));
            newresults(kk).rSquareStim(:,jj) = mean(results(kk).rSquareStim(:,INX{jj},2));
        end
    end
    % to do: add SE
end

% plot derived parameters
derivedTitles = {'time2peak', 'Rasymp'};
for kk = 1:nModels
    
    figure('Name', sprintf('%s %s', 'Derived parameters', func2str(results(kk).model))); hold on
    
    % plot explained variance
    subplot(1,3,1); plot(1:nChans,results(kk).rSquareConc, '.b', 'MarkerSize', 50, 'LineStyle', 'none')
    set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
    set(gca, 'Ylim', [0 1]);
    title('explained variance'); xlabel('visual area');  ylabel('R2'); set(gca, 'fontsize', 16);

    % plot derived parameters
    for p = 1:2
        subplot(1,3,p+1);
        plot(1:nChans,results(kk).derivedPrm(p,:), '.r', 'MarkerSize', 50, 'LineStyle', 'none')
        set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
        if p == 1, set(gca, 'Ylim', [0 0.5]), else, set(gca, 'YLim', [0 1]); end
        title(derivedTitles{p}); xlabel('visual area');  ylabel('parameter value'); set(gca, 'fontsize', 16);
    end 
    set(gcf, 'Position', [400 800 2000 600]);
end

% plot fitted parametesr
for kk = 1:nModels
    figure('Name', sprintf('%s %s', 'Fitted parameters', func2str(results(kk).model))); hold on
    
    % read in parameter files from json
    tmp = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', func2str(results(kk).model))));
    paramNames = strsplit(tmp.params,',');
    nParams = length(paramNames);
       
    % plot fitted parameters
    for p = 1:nParams
        subplot(2,ceil(nParams/2),p);
        plot(1:nChans,results(kk).params(p,:), '.k', 'MarkerSize', 50, 'LineStyle', 'none')
        set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
        title(paramNames{p}); xlabel('visual area'); ylabel('parameter value');set(gca, 'fontsize', 16);
    end
    set(gcf, 'Position', [400 200 2000 1200]);
end


end