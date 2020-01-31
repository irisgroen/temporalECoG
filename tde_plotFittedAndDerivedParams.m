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

% Currently assuming we'll never plot >6 models at a time
colors = {'r', 'b', 'c', 'm', 'g', 'y'}; 

% Extract timecourses and model names and derived parameters from results
modelNames = cell(1,nModels);
m = []; se = []; mp = [];

for kk = 1:nModels
    
    figure('Name', sprintf('Derived predictions %s', func2str(results(kk).model))); 
    
    for ii = 1:nChans
        subplot(ceil(sqrt(nChans)),ceil(sqrt(nChans)),ii); hold on
        l = cell(1,nModels);
        if ~dataWasAveraged
            %h = ciplot(m(kk,:,ii)-se(kk,:,ii), m(kk,:,ii)+se(kk,:,ii), [], colors{kk}, 0.25);
            %h = ciplot(se(kk,:,ii,1), se(kk,:,ii,2), [], colors{kk}, 0.25);
            %h.Annotation.LegendInformation.IconDisplayStyle = 'off';        
            %[m(kk,:,:), se(kk,:,:)] = averageAcrossAreas(results(kk).derivedPred, INX);
            plot(results(kk).derivedPred(:,INX{ii}), 'Color', colors{kk}, 'LineWidth', 2);        
            %mp = averageAcrossAreas(results(kk).derivedPrm, INX);
            %l{kk} = sprintf('%s median t2p = %0.2f median rasymp = %0.2f', ...
            %    func2str(results(kk).model), mp(kk,1,ii), mp(kk,2,ii));
        else
            plot(results(kk).derivedPred(:,ii), 'Color', colors{kk}, 'LineWidth', 2); 
            mp = results(kk).derivedPrm;
            l{kk} = sprintf('%s median t2p = %0.2f median rasymp = %0.2f', ...
                func2str(results(kk).model), mp(1,ii), mp(2,ii));
        end
        set(gca, 'Xlim', [0 1000]);
        modelNames{kk} = func2str(results(kk).model);
    end
    if dataWasAveraged, legend(l); end 
    title(channels.name{ii});
    set(gca, 'FontSize', 14);
    set(gcf, 'Position', [400 200 2000 1200]);

    % Save plot
    if saveFig
        figName = sprintf('derivedPredictions_%s', [modelNames{:}]);
        savePlot(figName, saveDir, dataWasAveraged)
    end
end

%% Plot derived parameters

% Separate figure for each model:
derivedTitles = {'time2peak', 'Rasymp'};

for kk = 1:nModels
    
    figure('Name', sprintf('Derived parameters %s', func2str(results(kk).model))); hold on
    
    % Plot explained variance
    subplot(1,3,1); hold on
    if ~dataWasAveraged
        [m, se, dat] = averageAcrossAreas(results(kk).rSquareConc, INX);
        for ii = 1:nChans, scatter(ones(1,size(dat{ii},2))*ii, dat{ii}, 40, [0.5 0.5 0.5], 'filled');end
        %errorbar(1:nChans, m, se, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
        errorbar(1:nChans, m, m-se(:,:,1), se(:,:,2)-m, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
    else
        m = results(kk).rSquareConc; 
        plot(1:nChans, m, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none')
    end

    set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
    set(gca, 'Ylim', [0 1]);
    title('explained variance'); xlabel('visual area');  ylabel('R2'); set(gca, 'fontsize', 16);

    % Plot derived parameters
    for p = 1:2
        subplot(1,3,p+1); hold on
        if ~dataWasAveraged
            [m, se, dat] = averageAcrossAreas(results(kk).derivedPrm(p,:), INX);
            for ii = 1:nChans, scatter(ones(1,size(dat{ii},2))*ii, dat{ii}, 40, [0.5 0.5 0.5], 'filled');end
            %errorbar(1:nChans, m, se, '.k', 'MarkerSize', 50, 'LineWidth', 2,'LineStyle', 'none', 'CapSize', 0)
            errorbar(1:nChans, m, m-se(:,:,1), se(:,:,2)-m, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
        else
            m = results(kk).derivedPrm(p,:); 
            plot(1:nChans, m, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none')
        end        
        set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
        if p == 1, set(gca, 'Ylim', [0 0.5]), else, set(gca, 'YLim', [0 1]); end
        title(derivedTitles{p}); xlabel('visual area');  ylabel('parameter value'); set(gca, 'fontsize', 16);
    end 
    set(gcf, 'Position', [400 800 2000 600]);
    
    % Save plot
    if saveFig
        figName = sprintf('derivedParams_%s', modelNames{kk});
        savePlot(figName, saveDir, dataWasAveraged)
    end
end

%% Plot fitted parameters

% Separate figure for each model:
for kk = 1:nModels
    
    figure('Name', sprintf('Fitted parameters %s', func2str(results(kk).model))); hold on
    
    % Read in parameter names from json
    tmp = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', func2str(results(kk).model))));
    paramNames = strsplit(tmp.params,',');
    nParams = length(paramNames);
       
    % Plot fitted parameters
    for p = 1:nParams
        subplot(2,ceil(nParams/2),p); hold on
        if ~dataWasAveraged
            [m, se, dat] = averageAcrossAreas(results(kk).params(p,:), INX);
            for ii = 1:nChans, scatter(ones(1,size(dat{ii},2))*ii, dat{ii}, 40, [0.5 0.5 0.5], 'filled');end
            %errorbar(1:nChans, m, se, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
            errorbar(1:nChans, m, m-se(:,:,1), se(:,:,2)-m, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
        else
            m = results(kk).params(p,:); 
            plot(1:nChans, m, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none')
        end
        set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
        title(paramNames{p}); xlabel('visual area'); ylabel('parameter value');set(gca, 'fontsize', 16);
    end
    set(gcf, 'Position', [400 200 2000 1200]);
    
    % Save plot
    if saveFig
        figName = sprintf('fittedParams_%s', modelNames{kk});
        savePlot(figName, saveDir, dataWasAveraged)
    end
end



end

% SUBROUTINES

% electrode averaging
function [m, se, indiv_points] = averageAcrossAreas(data, INX)
    nAreas = length(INX);  
    m = nan(size(data,1), nAreas); 
    se = nan(size(data,1), nAreas, 2);
    indiv_points = cell(nAreas,1);
    for jj = 1:nAreas
        elec_index = find(INX{jj});
        [m(:,jj),se(:,jj,:)] = calcmdsepct(data(:,elec_index),2);
        %m(:,jj)    = mean(data(:,elec_index),2);
        %se(:,jj)   = std(data(:,elec_index),0,2)/sqrt(length(elec_index));
        indiv_points{jj} = data(:,elec_index);
    end
end

% plot saving
function savePlot(figName, saveDir, dataWasAveraged)
    if ~dataWasAveraged
        figDir = fullfile(saveDir, 'individualelectrodes');
    else
        figDir = fullfile(saveDir, 'electrodeaverages');
    end
    if ~exist(figDir, 'dir'), mkdir(figDir), end
    saveas(gcf, fullfile(figDir, figName), 'png'); %close;
    %saveas(gcf, fullfile(figDir, figName), 'fig'); %close;
end
    
