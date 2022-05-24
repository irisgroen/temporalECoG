function [results] = tde_plotDerivedPredictions(results, channels, plotType, normalize, saveDir)

% Plot prediction for sustained stimulus
% plotType 1 = individual subplots per area, 2 = overlapping plots
% normalize (dividebymax), 1 = yes, 0 = no (default)
% 
% 2020 Iris Groen

if ~exist('plotType', 'var') || isempty(plotType), plotType = 2; end
if ~exist('normalize', 'var') || isempty(normalize), normalize = 0; end
if ~exist('saveDir', 'var') || isempty(saveDir), saveDir = []; end

%% Prep

nChans      = height(channels);
nModels     = size(results,2);

% Extract model names:
modelNames = cell(1,nModels);
for kk = 1:nModels
    modelNames{kk} = func2str(results(kk).model);
end

% Determine if data was averaged across elecs prior to fit; if not, average
% derived and fitted parameters now
if isfield(summary(channels), 'number_of_elecs')
    dataWasAveraged = true;
else
    dataWasAveraged = false;
    [chan_idx, channels] = groupElecsByVisualArea(channels);  
    nChans = height(channels);
end

% Are we saving figures?
if ~isempty(saveDir), saveFig = true; else, saveFig = false; end

%% Separate figure for each model and visual area:

if plotType == 1

    % Determine how many subplots to make
    nRow  = ceil(sqrt(nChans));
    nCol  = ceil(sqrt(nChans));
    if nChans <= (nRow*nCol)-nCol, nRow = nRow-1;end

    % Loop across models
    for kk = 1:nModels

        modelNames{kk} = func2str(results(kk).model);
        
        if normalize 
            figName = sprintf('DerivedPredictions normalized %s', modelNames{kk});
        else
            figName = sprintf('DerivedPredictions %s', modelNames{kk});
        end
            
        figure('Name', figName); 
        
        for ii = 1:nChans
            
            subplot(nRow,nCol,ii); hold on
            if ~dataWasAveraged
                %h = ciplot(m(kk,:,ii)-se(kk,:,ii), m(kk,:,ii)+se(kk,:,ii), [], colors{kk}, 0.25);
                %h = ciplot(se(kk,:,ii,1), se(kk,:,ii,2), [], colors{kk}, 0.25);
                %h.Annotation.LegendInformation.IconDisplayStyle = 'off';        
                %[m(kk,:,:), se(kk,:,:)] = averageAcrossAreas(results(kk).derivedPred, INX);
                m = results(kk).derived.pred_s(:,chan_idx(:,ii));
                %mp = averageAcrossAreas(results(kk).derivedPrm, INX);
                %l{kk} = sprintf('%s median t2p = %0.2f median rasymp = %0.2f', ...
                %    func2str(results(kk).model), mp(kk,1,ii), mp(kk,2,ii));
            else
                m = results(kk).derived.pred_s(:,ii);
                mp = results(kk).derived.params;
                l = sprintf('%s t2p = %0.2f rasymp = %0.2f', ...
                    func2str(results(kk).model), mp(1,ii), mp(2,ii));
            end
            if normalize, m = m/max(m); end
            plot(m, 'Color', 'k', 'LineWidth', 2); 
            set(gca, 'Xlim', [0 1000]);
            if normalize, set(gca, 'Ylim', [-0.1 1.1]);end

            title(channels.name{ii});
            if dataWasAveraged, legend(l); end 
            set(gca, 'FontSize', 14);
            if ii == 1, xlabel('Time'), ylabel('Predicted response');end
        end

        set(gcf, 'Position', [400 200 2000 1200]);

        % Save plot
        if saveFig
            figName = strrep(figName, ' ', '_');
            savePlot(figName, saveDir, dataWasAveraged)
        end
    end
end



    
%% One figure with subplots for each model, superimposed visual areas:

if plotType == 2
    
    if normalize 
        figName = sprintf('DerivedPredictions superimposed normalized %s', [modelNames{:}]);
    else
        figName = sprintf('DerivedPredictions superimposed %s', [modelNames{:}]);
    end
    colors = jet(nChans);
    
    figure('Name', figName); 

    % Determine how many subplots to make
    nRow  = ceil(sqrt(nModels));
    nCol  = ceil(sqrt(nModels));
    if nModels <= (nRow*nCol)-nCol, nRow = nRow-1; end

    % Loop across models
    for kk = 1:nModels
        modelNames{kk} = func2str(results(kk).model);
        subplot(nRow, nCol, kk); hold on;
        % Plot multiple areas together in one plot
        for ii = 1:nChans
            if ~dataWasAveraged
                m = mean(results(kk).derived.pred_s(:,chan_idx(:,ii)),2, 'omitnan');
                %se = std(results(kk).derived.pred(:,INX{ii}),0,2, 'omitnan');
                %h = ciplot(m-se, m+se, [], colors(ii,:), 0.25);
                %h.Annotation.LegendInformation.IconDisplayStyle = 'off';              
            else
                m = results(kk).derived.pred_s(:,ii);
                %m = squeeze(mean(results(kk).pred(:,[5 11],ii),2));
            end
            if normalize, m = m/max(m); end
            plot(m, 'Color', colors(ii,:), 'LineWidth', 2);        
        end
        set(gca, 'FontSize', 14);
        set(gca, 'Xlim', [0 1000]);
        if normalize, set(gca, 'Ylim', [-0.1 1.1]);end
        title(modelNames{kk});
        if kk ==1, legend(channels.name, 'Location', 'SouthEast'); xlabel('Time'), ylabel('Predicted response');end
    end

    set(gcf, 'Position', [400 200 2000 1200]);
    % Save plot
    if saveFig
        figName = strrep(figName, ' ', '_');
        savePlot(figName, saveDir, dataWasAveraged)
    end
end


end

% SUBROUTINES

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