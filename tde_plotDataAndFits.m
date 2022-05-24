function tde_plotDataAndFits(results, data, channels, stim_ts, stim_info, t, saveDir, conditionsOfInterest, timepointsOfInterest)

% Will generate plots of concatenated stimulus time courses (in gray),
% concatenated data time courses (in black) and model predictions (in colours).
%
% 2020 Iris Groen 

if ~exist('saveDir', 'var') || isempty(saveDir), saveDir = []; end
if ~exist('conditionsOfInterest', 'var') || isempty(conditionsOfInterest)
    conditionsOfInterest = {'CRF', 'ONEPULSE', 'TWOPULSE'};
end
if ~exist('timepointsOfInterest', 'var') || isempty(timepointsOfInterest)
    timepointsOfInterest = [t(1) t(end)];
end

nModels     = size(results,2);
nDatasets   = size(data,3);
nCond       = length(conditionsOfInterest);
%stim_info   = stim_info(contains(stim_info.name, conditionsOfInterest),:);
t_ind       = t>timepointsOfInterest(1) & t<=timepointsOfInterest(2);

% Determine if data was averaged across elecs prior to fit
if isfield(summary(channels), 'number_of_elecs')
    dataWasAveraged = true;
else
    dataWasAveraged = false;
end

%% Plot data and predictions
colors = {'r', 'b', 'c', 'm', 'g', 'y'}; % assuming we'll never plot >6 model fits at a time

% Prepare legend
l = cell(1,nModels+1);
l{1} = 'data';
for kk = 1:nModels, l{kk+1} = func2str(results(kk).model); end

if isfield(summary(channels), 'index_1')
	chan_index = channels.index_1;
else
    chan_index = channels.index;
end

% Loop over channels or channel averages
for ii = 1:nDatasets
    
    figure;
    
    d = data(t_ind,:,ii);
    %maxresp = max(d(:)); % scale stimulus to max across dataset
	maxresp = max(max(d(:,1:5))); % scale stimulus to max across CRF conditions

    % Loop over conditions 
    for jj = 1:length(conditionsOfInterest)
        subplot(nCond,1,jj); hold on
        
        inx = contains(stim_info.name, conditionsOfInterest{jj});
        cond = unique(stim_info.condition(inx));
        % plot stimulus
        h = plot(flatten(stim_ts(t_ind,inx))*maxresp, 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        
        % plot data
        plot(flatten(d(:,inx)), 'Color', 'k', 'LineWidth', 2); 
        titlestr = cell(1,nModels);
        
        % plot models
        for kk = 1:nModels           
            pred = results(kk).pred(t_ind,inx,ii);
            %plot(flatten(pred), 'Color', colors{kk}, 'LineStyle', '-.', 'LineWidth', 2);
            plot(flatten(pred), 'Color', colors{kk},  'LineWidth', 2);
            if isfield(results(kk).R2, 'concat_cond')
                R2val = mean(results(kk).R2.concat_cond(cond,ii));
            else
                R2val = mean(results(kk).R2.stim(inx,ii));
            end
            titlestr{kk} = sprintf('   r2 %s = %0.2f   ', func2str(results(kk).model), R2val);
        end
        
        % add title
        title(sprintf('%s: %s', conditionsOfInterest{jj}, [titlestr{:}]));
        
        % set axes
       % axis tight   
        set(gca, 'XTick',1:size(d,1):length(find(inx))*size(d,1), 'XTickLabel', []);
        axis tight
%         if contains(conditionsOfInterest{jj}, 'CRF')
%             set(gca, 'XTickLabel', stim_info.contrast(inx))
%         elseif contains(conditionsOfInterest{jj}, 'ONEPULSE')
%             set(gca, 'XTickLabel', stim_info.duration(inx))
%         else
%             set(gca, 'XTickLabel', stim_info.ISI(inx))
%         end
          %ylim([-0.2 1])
%         if ii == 1, ylim([-5 20]); end
%         if ii == 2, ylim([-5 15]); end
%         if ii == 3, ylim([-2 10]); end
%         if ii == 4, ylim([-0.5 5.5]); end
%         if ii == 5, ylim([-0.5 2]); end
%         if ii == 6, ylim([-0.5 2.5]); end

        set(gca, 'FontSize', 20);
        %xlabel('stimulus');
        %ylabel('broadband power');
        
        % add legend
        %if jj == 1, legend(l); end
    end
    set(gcf, 'Position', [400 200 1800 1200]);
    
    % Determine how to name the plot 
    if ~dataWasAveraged
        if isfield(summary(channels), 'benson14_varea') && isfield(summary(channels), 'wang15_mplbl')
            if isfield(summary(channels), 'subject_name')
                figureName = sprintf('%d_fits_%s_%s_%s_%s_%s', chan_index(ii), channels.benson14_varea{ii}, channels.wang15_mplbl{ii}, ...
                    channels.name{ii}, channels.subject_name{ii}, [l{2:end}]);
            else
                figureName = sprintf('%d_fits_%s_%s_%s_%s', chan_index(ii), channels.benson14_varea{ii}, channels.wang15_mplbl{ii}, ...
                    channels.name{ii}, [l{2:end}]);
            end
        else
            figureName = sprintf('fits_%s_%s', channels.name{ii}, [l{2:end}]);
        end
    else
        figureName = sprintf('fits_%s_%s', channels.name{ii}, [l{2:end}]);
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