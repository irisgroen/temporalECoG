function [results] = tde_plotParams(results, channels, saveDir)

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

%% Plot derived parameters

% Separate figure for each model:
derivedTitles = [{'R2 concatenated'} results(1).derived.names];
m_all = [];
se_all = [];

for kk = 1:nModels
    
    figure('Name', sprintf('Derived parameters %s', modelNames{kk})); hold on
    set(gcf, 'Position', [400 200 2000 1200]);
    
    % Plot explained variance
    subplot(2,3,1); hold on
    if ~dataWasAveraged
        [m, se, dat] = averageAcrossAreas(results(kk).R2.concat_all, INX);
        for ii = 1:nChans, scatter(ones(1,size(dat{ii},2))*ii, dat{ii}, 40, [0.5 0.5 0.5], 'filled');end
        %errorbar(1:nChans, m, se, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
        errorbar(1:nChans, m, m-se(:,:,1), se(:,:,2)-m, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
        se_all(1,:,:,kk) = squeeze(se);
    else
        m = results(kk).R2.concat_all; 
        plot(1:nChans, m, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none')
    end
    m_all(1,:,kk) = squeeze(m);
    
    set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
    set(gca, 'Ylim', [0 1]);
    title(derivedTitles{1}); xlabel('visual area');  ylabel('R2'); set(gca, 'fontsize', 16);

    % Plot derived parameters
    subplotinx = [2 3 5 6];
    for p = 1:4
        subplot(2,3,subplotinx(p)); hold on
        if ~dataWasAveraged
            [m, se, dat] = averageAcrossAreas(results(kk).derived.params{p}, INX);
            for ii = 1:nChans, scatter(ones(1,size(dat{ii},2))*ii, dat{ii}, 40, [0.5 0.5 0.5], 'filled');end
            %errorbar(1:nChans, m, se, '.k', 'MarkerSize', 50, 'LineWidth', 2,'LineStyle', 'none', 'CapSize', 0)
            errorbar(1:nChans, m, m-se(:,:,1), se(:,:,2)-m, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
            se_all(p+1,:,:,kk) = squeeze(se);
        else
            m = results(kk).derived.params{p};
            if size(m,2) > 1
                cmap = num2cell(flipud(gray(size(m,2)+1)),2);
                h = plot(1:nChans, m, 'k.-', 'MarkerSize', 20, 'LineWidth', 2);
                set(h, {'color'}, cmap(2:end));
                %if p == 4, legend({'80%', '90%', '95%', '100%'},'Location', 'NorthWest'); end
                if p == 4, m = m(:,3); end % 95%
                if p == 3, m = m(:,4); end % 0.267s vs 0.133s 
            else
                plot(1:nChans, m, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none')
            end
        end        
        set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
        if p == 1, set(gca, 'Ylim',[0 0.2]); end
        if p == 2, set(gca, 'Ylim',[0 1]); end
        if p == 3, set(gca, 'Ylim',[0.5 1.2]); end
        if p == 4, set(gca, 'Ylim',[0 1.2]); end
        title(derivedTitles{p+1}); xlabel('visual area');  ylabel('parameter value'); set(gca, 'fontsize', 16);
        m_all(p+1,:,kk) = squeeze(m);
    end
    
    % Save plot
    if saveFig
        figName = sprintf('derivedParams_%s', modelNames{kk});
        savePlot(figName, saveDir, dataWasAveraged)
    end
end

% All models together in one plot:

figure('Name', 'Derived parameters - all models'); hold on
set(gcf, 'Position', [400 800 2000 600]);

%colors = parula(nModels);

for jj = 1:size(m_all,1)

    subplot(1,5,jj); hold on
    m = squeeze(m_all(jj,:,:));
    
    if size(m,2) == 1 % If there's just one area, add a dummy column to make sure bars will still be grouped
        m=cat(2,m,nan(size(m)))';
    end
    h = bar(m);
    set(h, 'BarWidth', 1); 

    numgroups = size(m,1);
    numbars = size(m,2);
    groupwidth = min(0.8,numbars/(numbars+1.5));

    if ~isempty(se_all)
        se = squeeze(se_all(jj,:,:,:));
        if nModels > 1
            neg = m-se(:,1,:);
            pos = se(:,2,:)-m;
            for ii = 1:numbars
                x = (1:numgroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar    
                errorbar(x, m(:,ii), neg(:,ii), pos(:,ii), 'k', 'LineWidth', 2,  'LineStyle', 'none', 'CapSize', 0);
            end
        else
            neg = m - se(:,1)';
            pos = se(:,2)'-m;
            x = 1:length(m);
            errorbar(x, m, neg, pos, 'k', 'LineWidth', 2,  'LineStyle', 'none', 'CapSize', 0);
        end 
    end
    set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
    title(derivedTitles{jj}); xlabel('visual area'); 
    if jj == 1, legend(modelNames); set(gca, 'Ylim', [0 1]); end
    if jj == 2, set(gca, 'Ylim',[0 0.2]); end
    if jj == 3, set(gca, 'Ylim',[0 1]); end
    if jj == 4, set(gca, 'Ylim',[0.5 1]); end
	if jj == 5, set(gca, 'Ylim',[0 1]); end
    set(gca, 'fontsize', 16);
end

% Save plot
if saveFig
    figName = sprintf('derivedParams_%s_allmodels', [modelNames{:}]);
    savePlot(figName, saveDir, dataWasAveraged)
end

%% Plot fitted parameters

% Separate figure for each model:
for kk = 1:nModels
    
    figure('Name', sprintf('Fitted parameters %s', modelNames{kk})); hold on
    
    % Read in parameter names from json
    tmp = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', modelNames{kk})));
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
    
