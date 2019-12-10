function [results] = tde_plotFittedAndDerivedParams(results, channels, saveDir)

if ~exist('saveDir', 'var') || isempty(saveDir), saveDir = []; end

nChans      = height(channels);
nModels     = size(results,2);

% Determine if data was averaged across elecs prior to fit; if not, average
% derived and fitted parameters now
if isfield(summary(channels), 'number_of_elecs')
    dataWasAveraged = true;
else
    dataWasAveraged = false;
    [INX, channels] = groupElecsByVisualArea(channels);  
    nChans = size(INX,2);
end

% Determine if we're saving figures
if ~isempty(saveDir), saveFig = true; else, saveFig = false; end

%% Plot prediction for sustained stimulus

% Plot multiple models in each subplot

% Currently assuming we'll never plot >4 models at a time
colors = {'r', 'b', 'c', 'm'}; 

% Extract timecourses and model names and derived parameters from results
modelNames = [];
m = []; se = []; mp = [];
for kk = 1:nModels
    if ~dataWasAveraged
        [m(kk,:,:), se(kk,:,:)] = averageAcrossAreas(results(kk).derivedPred, INX);
        [mp(kk,:,:)] = averageAcrossAreas(results(kk).derivedPrm, INX);
    else
        m(kk,:,:) = results(kk).derivedPred; se = [];
    end
    modelNames{kk} = func2str(results(kk).model);
end

% Make plot
figure('Name', sprintf('%s', 'Derived predictions')); 
for ii = 1:nChans
    subplot(ceil(sqrt(nChans)),ceil(sqrt(nChans)),ii); hold on
    l = cell(1,nModels);
    for kk = 1:nModels
        if exist('se', 'var')
            h = ciplot(m(kk,:,ii)-se(kk,:,ii), m(kk,:,ii)+se(kk,:,ii), [], colors{kk}, 0.25);
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';        
        end
        plot(m(kk,:,ii), 'Color', colors{kk}, 'LineWidth', 2);        
        set(gca, 'Xlim', [0 1000]);
        l{kk} = sprintf('%s t2p = %0.2f rasymp = %0.2f', ...
            func2str(results(kk).model), mp(kk,1,ii), mp(kk,2,ii));
    end
	legend(l); 
    title(channels.name{ii});
    set(gca, 'FontSize', 14);
    
end
set(gcf, 'Position', [400 200 2000 1200]);

% Save plot
if saveFig
    figName = sprintf('derivedPredictions_model%s', [modelNames{:}]);
    savePlot(figName, saveDir, dataWasAveraged)
end

%% Plot derived parameters

% Separate figure for each model:
derivedTitles = {'time2peak', 'Rasymp'};
for kk = 1:nModels
    
    figure('Name', sprintf('%s %s', 'Derived parameters', func2str(results(kk).model))); hold on
    
    % Plot explained variance
    subplot(1,3,1); 
    if ~dataWasAveraged
        [m, se] = averageAcrossAreas(results(kk).rSquareConc, INX);
    else
        m = results(kk).rSquareConc; se = [];
    end
    errorbar(1:nChans, m, se, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
    set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
    set(gca, 'Ylim', [0 1]);
    title('explained variance'); xlabel('visual area');  ylabel('R2'); set(gca, 'fontsize', 16);

    % Plot derived parameters
    for p = 1:2
        subplot(1,3,p+1);
        if ~dataWasAveraged
            [m, se] = averageAcrossAreas(results(kk).derivedPrm(p,:), INX);
        else
            m = derivedPrm(p,:); se = [];
        end
        errorbar(1:nChans, m, se, '.k', 'MarkerSize', 50, 'LineWidth', 2,'LineStyle', 'none', 'CapSize', 0)
        set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
        if p == 1, set(gca, 'Ylim', [0 0.5]), else, set(gca, 'YLim', [0 1]); end
        title(derivedTitles{p}); xlabel('visual area');  ylabel('parameter value'); set(gca, 'fontsize', 16);
    end 
    set(gcf, 'Position', [400 800 2000 600]);
    
    % Save plot
    if saveFig
        figName = sprintf('derivedParams_model%s', [modelNames{:}]);
        savePlot(figName, saveDir, dataWasAveraged)
    end
end


%% Plot fitted parameters

% Separate figure for each model:
for kk = 1:nModels
    figure('Name', sprintf('%s %s', 'Fitted parameters', func2str(results(kk).model))); hold on
    
    % Read in parameter names from json
    tmp = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', func2str(results(kk).model))));
    paramNames = strsplit(tmp.params,',');
    nParams = length(paramNames);
       
    % Plot fitted parameters
    for p = 1:nParams
        subplot(2,ceil(nParams/2),p);
        if ~dataWasAveraged
            [m, se] = averageAcrossAreas(results(kk).params(p,:), INX);
        else
            m = results(kk).params(p,:); se = [];
        end
        errorbar(1:nChans, m, se, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
        set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
        title(paramNames{p}); xlabel('visual area'); ylabel('parameter value');set(gca, 'fontsize', 16);
    end
    set(gcf, 'Position', [400 200 2000 1200]);
    
    % Save plot
    if saveFig
        figName = sprintf('fittedParams_model%s', [modelNames{:}]);
        savePlot(figName, saveDir, dataWasAveraged)
    end
end



end

% SUBROUTINES

% electrode averaging
function [m, se] = averageAcrossAreas(data, INX)
    nAreas = length(INX);  
    m = nan(size(data,1), nAreas); 
    se = nan(size(data,1), nAreas);
    for jj = 1:nAreas
        elec_index = find(INX{jj});
        m(:,jj)    = mean(data(:,elec_index),2);
        se(:,jj)   = std(data(:,elec_index),0,2)/sqrt(length(elec_index));
    end
end

% plot saving
function savePlot(figName, saveDir, dataWasAveraged)
    if ~dataWasAveraged
        figDir = fullfile(saveDir, 'individualelectrodes');
    else
        figDir = fullfile(saveDir, 'electrodeaverages');
    end
    saveas(gcf, fullfile(figDir, figName), 'png'); close;
end
    
