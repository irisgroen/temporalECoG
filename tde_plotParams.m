function [results] = tde_plotParams(results, channels, saveDir, opts)

% Plot derived and fitted parameters for models in results
%
% 2020 Iris Groen

if ~exist('saveDir', 'var') || isempty(saveDir), saveDir = []; end
if ~exist('opts', 'var') || isempty(opts), opts = struct(); end

if ~isfield(opts, 'plotindivpoints'), opts.plotindivpoints = false;end

%% Prep

nChans      = height(channels);
nModels     = size(results,2);

% Extract model names:
modelNames = cell(1,nModels);
for kk = 1:nModels
    modelNames{kk} = func2str(results(kk).model);
end

% Determine if data was averaged across elecs prior to fit
if isfield(summary(channels), 'number_of_elecs')
    dataWasAveraged = true;
else
    dataWasAveraged = false;
    %[chan_idx1, channels1] = groupElecsByVisualArea(channels, 'fixedassignment');  
	%[chan_idx2, channels2] = groupElecsByVisualArea(channels, 'probabilisticresample');  
    %channels = channels1;
    [~, channels, group_prob] = groupElecsByVisualArea(channels, 'probabilisticresample');
    nChans = height(channels);
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
    
    % Collect the data to plot
    data{1} = results(kk).R2.concat_all;
    for p = 1:length(derivedTitles)-1
       data{p+1} = results(kk).derived.params(p,:);
    end
%     for p = 1:3
%         data{p} = results(kk).R2.concat_cond(p,:);
%     end
       
    % Plot 
    subplotinx = 1:length(derivedTitles);%^[1 2 3 5 6];
    for dd = 1:length(data)
        subplot(1,length(derivedTitles),subplotinx(dd)); hold on
        if ~dataWasAveraged
            [m, se] = averageWithinArea(data{dd}, group_prob);
            if opts.plotindivpoints
                for ii = 1:nChans, scatter(ones(1,size(dat{ii},2))*ii, dat{ii}, 30, [0.7 0.7 0.7], 'filled');end
            end
            errorbar(1:nChans, m, m-se(:,1)', se(:,2)'-m, '.k', 'MarkerSize', 30, 'LineWidth', 4, 'LineStyle', 'none', 'CapSize', 0)
            %errorbar(1:nChans, m, m-se(:,:,1), se(:,:,2)-m, '.b', 'MarkerSize', 30, 'LineWidth', 4, 'LineStyle', 'none', 'CapSize', 0)
            %[m, se] = averageAcrossElecsWithinArea(data{dd}, chan_idx2);
            %errorbar(1:nChans, m, m-se(:,:,1), se(:,:,2)-m, '.c', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
            se_all(dd,:,:,kk) = squeeze(se);
        else
            %m = results(kk).R2.concat_all; 
            m = data{dd};
            if size(m,1) > 1
                cmap = num2cell(flipud(gray(size(m,1)+1)),2);
                h = plot(1:nChans, m, 'k.-', 'MarkerSize', 20, 'LineWidth', 2);
                set(h, {'color'}, cmap(2:end));
                %if p == 4, legend({'80%', '90%', '95%', '100%'},'Location', 'NorthWest'); end
                if p == 4, m = m(3,:); end % 95%
                if p == 3, m = m(4,:); end % 0.267s vs 0.133s 
            else
                plot(1:nChans, m, '.k', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none')
            end
        end
        m_all(dd,:,kk) = squeeze(m);
        
        % set axes and axes labels
        set(gca, 'Xlim', [0 nChans+1], 'XTick', 1:nChans, 'XTickLabel', channels.name, 'XTickLabelRotation', 45);
        set(gca, 'Ylim', [0 1]);
        title(derivedTitles{dd}); xlabel('visual area'); set(gca, 'fontsize', 16);
        if dd == 1, set(gca, 'Ylim',[0 1]); end
        if dd == 2, set(gca, 'Ylim',[0 0.2]); end
        if dd == 3, set(gca, 'Ylim',[0 1]); end
        if dd == 4, set(gca, 'Ylim',[0.5 1.5]); end
        if dd == 5, set(gca, 'Ylim',[0 1.2]); end
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
nSubPlot = 3;%size(m_all,1)
for jj = 1:nSubPlot

    subplot(1,nSubPlot,jj); hold on
    %figure('Name', sprintf('allModels %s', derivedTitles{jj}));hold on;
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
            neg = m-squeeze(se(:,1,:));
            pos = squeeze(se(:,2,:))-m;
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
    if jj == 1, set(gca, 'Ylim', [0 1]); end
    if jj == 2, set(gca, 'Ylim',[0 1]); end
    if jj == 3, set(gca, 'Ylim',[0 1]); end
    if jj == 4, set(gca, 'Ylim',[0.5 1]); end
	if jj == 5, set(gca, 'Ylim',[0 1]); end
    if jj == nSubPlot, legend(modelNames, 'Location', 'NorthEast'); end
    set(gca, 'fontsize', 16);
    
    % Save plot
    if saveFig
        figName = strrep(get(gcf,'Name'),' ','_');
        savePlot(figName, saveDir, dataWasAveraged)
    end
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
            %[m, se, dat] = averageAcrossElecsWithinArea(results(kk).params(p,:), chan_idx1);
            [m, se] = averageWithinArea(results(kk).params(p,:), group_prob);
            if opts.plotindivpoints
                for ii = 1:nChans, scatter(ones(1,size(dat{ii},2))*ii, dat{ii}, 30, [0.7 0.7 0.7], 'filled');end
            end
            errorbar(1:nChans, m, m-se(:,1)', se(:,2)'-m, '.b', 'MarkerSize', 30, 'LineWidth', 4, 'LineStyle', 'none', 'CapSize', 0)
            %errorbar(1:nChans, m, se, '.k', 'MarkerSize', 50, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
            %errorbar(1:nChans, m, m-se(:,:,1), se(:,:,2)-m, '.b', 'MarkerSize', 30, 'LineWidth', 4, 'LineStyle', 'none', 'CapSize', 0)
            %[m, se] = averageAcrossElecsWithinArea(results(kk).params(p,:), chan_idx2);
            %errorbar(1:nChans, m, m-se(:,:,1), se(:,:,2)-m, '.c', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
        else
            m = results(kk).params(p,:); 
            plot(1:nChans, m, '.k', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none')
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
    if ~exist(saveDir), mkdir(saveDir); end
    if ~dataWasAveraged
        figDir = fullfile(saveDir, 'individualelectrodes');
    else
        figDir = fullfile(saveDir, 'electrodeaverages');
    end
    if ~exist(figDir, 'dir'), mkdir(figDir), end
    saveas(gcf, fullfile(figDir, figName), 'png'); %close;
    %saveas(gcf, fullfile(figDir, figName), 'fig'); %close;
end
    
