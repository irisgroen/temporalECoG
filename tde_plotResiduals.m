function tde_plotResiduals(results, data, channels, stim_ts, stim_info, t, conditionsOfInterest, saveDir)

% Plot concatenated residual time courses for models in results
%
% 2020 Iris Groen

if ~exist('conditionsOfInterest', 'var') || isempty(conditionsOfInterest)
    conditionsOfInterest = {'CRF', 'ONEPULSE', 'TWOPULSE'};
end
if ~exist('saveDir', 'var') || isempty(saveDir), saveDir = []; end

nModels     = size(results,2);
nDatasets   = size(data,3);
nCond       = length(conditionsOfInterest);
stim_info   = stim_info(contains(stim_info.name, conditionsOfInterest),:);

% Determine if data was averaged across elecs prior to fit; if not, average
% derived and fitted parameters now
if isfield(summary(channels), 'number_of_elecs')
    dataWasAveraged = true;
else
    dataWasAveraged = false;
end

%% plot data and predictions
colors = {'r', 'b', 'c', 'm', 'g', 'y'}; % assuming we'll never plot >6 model fits at a time

% Prepare legend
l = cell(1,nModels);
for kk = 1:nModels, l{kk} = func2str(results(kk).model); end

% Loop over channels or channel averages
for ii = 1:3%nDatasets
    
    figure;
    
    d = data(:,:,ii);
    maxresp = max(d(:)); % scale stimulus to max across dataset
    
    % Loop over conditions 
    for jj = 1:length(conditionsOfInterest)
        subplot(nCond,1,jj); hold on
        inx = contains(stim_info.name, conditionsOfInterest{jj});
        
        % plot stimulus
        h = plot(flatten(stim_ts(:,inx))*maxresp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        % plot residuals
        for kk = 1:nModels           
            pred = results(kk).pred(:,inx,ii);
            res = d(:,inx)-pred;
            plot(flatten(res), 'Color', colors{kk}, 'LineWidth', 2);
            if isfield(results(kk).R2, 'concat_cond')
                R2val = mean(results(kk).R2.concat_cond(jj,ii));
            else
                R2val = mean(results(kk).R2.stim(inx,ii));
            end
            titlestr{kk} = sprintf('   r2 %s = %0.2f   ', func2str(results(kk).model), R2val);
        end
        
        % add title
        title(sprintf('%s: %s', conditionsOfInterest{jj}, [titlestr{:}]));
        
        % set axes
        axis tight   
        set(gca, 'XTick',1:length(t):length(find(inx))*length(t), 'XTickLabel', []);
%         if contains(conditionsOfInterest{jj}, 'CRF')
%             set(gca, 'XTickLabel', stim_info.contrast(inx))
%         elseif contains(conditionsOfInterest{jj}, 'ONEPULSE')
%             set(gca, 'XTickLabel', stim_info.duration(inx))
%         else
%             set(gca, 'XTickLabel', stim_info.ISI(inx))
%         end
        if ii == 1, ylim([-5 25]); end
        if ii == 2, ylim([-5 15]); end
        if ii == 3, ylim([-2 8]); end

        set(gca, 'FontSize', 14);
        xlabel('stimulus');
        ylabel('response');
        
        % add legend
        if jj == 1, legend(l); end
    end
    set(gcf, 'Position', [400 200 1800 1200]);
    
    % Determine how to name the plot 
    if ~dataWasAveraged
        if isfield(summary(channels), 'bensonarea') 
            if isfield(summary(channels), 'subject_name')
                figureName = sprintf('residuals_%s_%s_%s_%s_%s', channels.bensonarea{ii}, channels.wangarea{ii}, ...
                    channels.name{ii}, channels.subject_name{ii}, [l{:}]);
            else
                figureName = sprintf('residuals_%s_%s_%s_%s', channels.bensonarea{ii}, channels.wangarea{ii}, ...
                    channels.name{ii}, [l{:}]);
            end
        else
            figureName = sprintf('residuals_%s_%s', channels.name{ii}, [l{:}]);
        end
    else
        figureName = sprintf('residuals_%s_%s', channels.name{ii}, [l{:}]);
    end
    
    set(gcf, 'Name', figureName);
    
    % Determine whether to save it and if so where
    if ~isempty(saveDir)
        if ~dataWasAveraged        
            figDir = fullfile(saveDir, 'individualelectrodes');
        else
            figDir = fullfile(saveDir, 'electrodeaverages');
        end
        if ~exist(figDir, 'dir'), mkdir(figDir), end
        saveas(gcf, fullfile(figDir, figureName), 'png'); 
        %saveas(gcf, fullfile(figDir, figureName), 'fig'); 
        close;
    end
end

end