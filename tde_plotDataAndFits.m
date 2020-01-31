function tde_plotDataAndFits(results, data, channels, stim_ts, stim_info, t, conditionsOfInterest, saveDir)

% potentially change this to tde_generateReport 

if ~exist('conditionsOfInterest', 'var') || isempty(conditionsOfInterest)
    conditionsOfInterest = {'CRF', 'ONEPULSE', 'TWOPULSE'};
end
if ~exist('saveDir', 'var') || isempty(saveDir), saveDir = []; end

nModels     = size(results,2);
nDatasets   = size(data,3);
nCond       = length(conditionsOfInterest);
stim_info   = stim_info(contains(stim_info.name, conditionsOfInterest),:);
%% plot data and predictions
colors = {'r', 'b', 'c', 'm', 'g', 'y'}; % assuming we'll never plot >6 model fits at a time

% Prepare legend
l = cell(1,nModels+1);
l{1} = 'data';
for kk = 1:nModels, l{kk+1} = func2str(results(kk).model); end

% Loop over channels or channel averages
for ii = 1:nDatasets
    
    figure('Name', sprintf('%s %s', 'Predictions', channels.name{ii}));
    
    d = data(:,:,ii);
    maxresp = max(d(:)); % scale stimulus to max across dataset
    
    % Loop over conditions 
    for jj = 1:length(conditionsOfInterest)
        subplot(nCond,1,jj); hold on
        inx = contains(stim_info.name, conditionsOfInterest{jj});
        % plot stimulus
        h = plot(flatten(stim_ts(:,inx))*maxresp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        % plot data
        plot(flatten(d(:,inx)), 'Color', 'k', 'LineWidth', 2); 
        titlestr = cell(1,nModels);
        
        % plot models
        for kk = 1:nModels           
            pred = results(kk).pred(:,inx,ii);
            plot(flatten(pred), 'Color', colors{kk}, 'LineWidth', 2);
            titlestr{kk} = sprintf('   r2 %s = %0.2f   ', func2str(results(kk).model), mean(results(kk).rSquareStim(inx,ii)));
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
        set(gca, 'FontSize', 16);
        xlabel('stimulus');
        ylabel('response');
        
        % add legend
        if jj == 1, legend(l); end
    end
    set(gcf, 'Position', [400 200 1800 1200]);

    if ~isempty(saveDir)
        if ~isfield(summary(channels), 'number_of_elecs') % this is non-averaged data
            if isfield(summary(channels), 'bensonarea') 
                if isfield(summary(channels), 'subject_name')
                    figureName = sprintf('%s_%s_%s_%s_%s', channels.bensonarea{ii}, channels.wangarea{ii}, ...
                        channels.name{ii}, channels.subject_name{ii}, [l{2:end}]);
                else
                    figureName = sprintf('%s_%s_%s_%s', channels.bensonarea{ii}, channels.wangarea{ii}, ...
                        channels.name{ii}, [l{2:end}]);
                end
            else
                figureName = sprintf('%s_%s', ...
                    channels.name{ii}, [l{2:end}]);
            end
            figDir = fullfile(saveDir, 'individualelectrodes');
        else
            figureName = sprintf('%s_%s', channels.name{ii}, [l{2:end}]);
            figDir = fullfile(saveDir, 'electrodeaverages');
        end
        if ~exist(figDir, 'dir'), mkdir(figDir), end
        saveas(gcf, fullfile(figDir, figureName), 'png'); 
        %saveas(gcf, fullfile(figDir, figureName), 'fig'); 
        close;
    end
    close;
end

end