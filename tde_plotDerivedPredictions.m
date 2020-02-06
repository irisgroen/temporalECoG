function [results] = tde_plotDerivedPredictions(results, channels, saveDir)

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
    [INX, channels] = groupElecsByVisualArea(channels);  
    nChans = size(INX,2);
end

% Are we saving figures?
if ~isempty(saveDir), saveFig = true; else, saveFig = false; end

%% Plot prediction for sustained stimulus

% Separate figure for each model:

% Extract timecourses and model names and derived parameters from results

for kk = 1:nModels
    
	modelNames{kk} = func2str(results(kk).model);
    figure('Name', sprintf('Derived predictions %s', modelNames{kk})); 
    
    for ii = 1:nChans
        subplot(ceil(sqrt(nChans)),ceil(sqrt(nChans)),ii); hold on
        l = cell(1,nModels);
        if ~dataWasAveraged
            %h = ciplot(m(kk,:,ii)-se(kk,:,ii), m(kk,:,ii)+se(kk,:,ii), [], colors{kk}, 0.25);
            %h = ciplot(se(kk,:,ii,1), se(kk,:,ii,2), [], colors{kk}, 0.25);
            %h.Annotation.LegendInformation.IconDisplayStyle = 'off';        
            %[m(kk,:,:), se(kk,:,:)] = averageAcrossAreas(results(kk).derivedPred, INX);
            plot(results(kk).derived.pred(:,INX{ii}), 'Color', 'k', 'LineWidth', 2);        
            %mp = averageAcrossAreas(results(kk).derivedPrm, INX);
            %l{kk} = sprintf('%s median t2p = %0.2f median rasymp = %0.2f', ...
            %    func2str(results(kk).model), mp(kk,1,ii), mp(kk,2,ii));
        else
            plot(results(kk).derived.pred(:,ii), 'Color', 'k', 'LineWidth', 2); 
            mp = results(kk).derived.params;
            l{kk} = sprintf('%s median t2p = %0.2f median rasymp = %0.2f', ...
                func2str(results(kk).model), mp(1,ii), mp(2,ii));
        end
        set(gca, 'Xlim', [0 1000]);
        
        title(channels.name{ii});
    end
    
    %if dataWasAveraged, legend(l); end 
    set(gca, 'FontSize', 14);
    set(gcf, 'Position', [400 200 2000 1200]);

    % Save plot
    if saveFig
        figName = sprintf('derivedPredictions_%s', [modelNames{kk}]);
        savePlot(figName, saveDir, dataWasAveraged)
    end
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