function tde_fitPRFs(data, bar_apertures, opt, doPlots, saveDir, resultsStr)

% <saveDir> path to save parameters and fits
%   default: fullfile(analysisRootPath, 'prfs')
% <resultsStr> string to add to the save filename for the results.mat
%   created by analyzePRF, if results are saved
%   default: 'prfs'

% <doPlots>
if ~exist('doPlots','var') || isempty(doPlots)
    doPlots = false; % boolean
end

% <saveDir>
if ~exist('saveDir', 'var'), saveDir = fullfile(analysisRootPath, 'prfs'); end

% <resultsStr>
if ~exist('resultsStr', 'var'), resultsStr = 'prfs'; end


if doPlots
    plotSaveDir = fullfile(analysisRootPath, 'figures', 'prfs');
    if ~exist(plotSaveDir, 'dir'); mkdir(fullfile(plotSaveDir));end
end

tr = 1; % no HRF for ECoG data

% Fit PRF models for each electrode in each subject
nSubjects = length(data);

% Loop over subjects
for ii = 1:nSubjects

    subject = data{ii}.subject;
	channels = data{ii}.channels;
        
    if ~isempty(data)
        
        % Average runs
        data2fit = mean(data{ii}.ts,3);
    
        % Define stimulus
        stim_inx = data{ii}.stim_inx;
        stimulus = bar_apertures(:,:,stim_inx);
        %stimulus = {bar_apertures,bar_apertures,bar_apertures,bar_apertures};

        results = analyzePRF(stimulus,data2fit, tr, opt);

%         % Save fits to results directory
%         if ~isempty(saveDir)
% 
%             if ~exist(saveDir, 'dir'); mkdir(saveDir); end
% 
%             saveName = sprintf('%s_%s', subject, resultsStr);
%             saveName = fullfile(saveDir, saveName);
%             fprintf('[%s] Saving results to %s \n', mfilename, saveName);
% 
%             if exist(sprintf('%s.mat',saveName),'file')
%                 warning('[%s] Results file already exists! Writing new file with date-time stamp.',mfilename);
%                 saveName = sprintf('%s_%s', saveName, datestr(now,30));
%                 fprintf('[%s] Saving results to %s \n', mfilename, saveName);
%             end
%             save(saveName, 'channels','data','results','tr','opt','subject');  
%         end
    
        % Make plots of the estimated PRFs and PRF fits

        if doPlots
            
            channels = data{ii}.channels;
            nChans = height(channels);
            figureName = sprintf('%s_prftimeseriesfits', subject);

            % Timeseries + fits
            ecog_plotPRFtsfits(data2fit, stimulus, results, channels)
            set(gcf, 'Name', figureName);
            set(gcf, 'Position', get(0,'screensize'));
            set(gca, 'FontSize', 14);
            saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;

            % PRFs
            coloropt = 0;
            figureName = sprintf('%s_prfs', subject);
            
            ecog_plotPRFs(results, stimulus, channels, coloropt)  
            set(gcf, 'Position', get(0,'screensize'));
            set(gca, 'FontSize', 14);
            saveas(gcf, fullfile(plotSaveDir, figureName), 'png'); close;
            
        end
    end
end

end