function [out] = tde_fitPRFs(data, bar_apertures, opt, doPlots, saveDir, resultsStr, plotSaveDir)

% Fits a pRF model to PRF time series using a modified version of Kendrick
% Kay's analyzePRF.m, created by Ken Yuasa (see
% https://github.com/WinawerLab/ECoG_utils/)
% 
% <saveDir> path to save parameters and fits
%   default: fullfile(analysisRootPath, 'prfs')
% <resultsStr> string to add to the save filename for the results.mat
%   created by analyzePRF, if results are saved
%   default: 'prfs'
% 
% 2022 Iris Groen

% <doPlots>
if ~exist('doPlots','var') || isempty(doPlots), doPlots = false; end

% <saveDir>
if ~exist('saveDir', 'var') || isempty(saveDir), saveDir = fullfile(analysisRootPath, 'prfs'); end

% <resultsStr>
if ~exist('resultsStr', 'var') || isempty(resultsStr), resultsStr = 'prffits'; end

% <plotSaveDir>
if ~exist('plotSaveDir','var') || isempty(plotSaveDir)
    plotSaveDir = fullfile(analysisRootPath, 'figures', 'prfs');
end
if ~exist(plotSaveDir, 'dir'); mkdir(fullfile(plotSaveDir));end

% Resize images to speed up calculations
bar_apertures = imresize(bar_apertures, [100 100], 'nearest');

tr = 1; % no HRF for ECoG data

% Fit PRF models for each electrode in each subject
nSubjects = length(data);

out = cell(nSubjects,1);

% Loop over subjects
for ii = 1:nSubjects

    subject = data{ii}.subject;
	channels = data{ii}.channels;
    %data2fit = [];
    %stimulus = [];
    
    if ~isempty(data)
        
        % Average runs
        if ~isfield(data{ii}, 'ts')
            continue
        else
            
            % Define stimulus and data
%            stim_inx = data{ii}.stim_inx;
%             for jj = 1:size(data{ii}.ts,3)
%                 data2fit{jj} = data{ii}.ts(:,:,jj);
%                 stimulus{jj} = double(bar_apertures(:,:,stim_inx));
%             end
%             
            data2fit = mean(data{ii}.ts,3);
            stimulus = {bar_apertures(:,:,data{ii}.stim_inx)};
            
            %results = analyzePRF_bounds(stimulus, data2fit, tr, opt);
            results = analyzePRFdog(stimulus, data2fit, tr, opt);
            results.channels = channels;
            results.subject  = subject;
            out{ii}          = results;
            
            % Save fits to results directory
            if ~isempty(saveDir)
    
                if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    
                saveName = sprintf('sub-%s_%s', subject, resultsStr);
                saveName = fullfile(saveDir, saveName);
                fprintf('[%s] Saving results to %s \n', mfilename, saveName);
    
                if exist(sprintf('%s.mat',saveName),'file')
                    warning('[%s] Results file already exists! Writing new file with date-time stamp.',mfilename);
                    saveName = sprintf('%s_%s', saveName, datestr(now,30));
                    fprintf('[%s] Saving results to %s \n', mfilename, saveName);
                end
                save(saveName, 'data2fit','results', 'stimulus');  
            end

            % Make plots of the estimated PRFs and PRF fits

            if doPlots
                close all;
                channels = data{ii}.channels;
                f_ind = checkForHDgrid(channels);

                % Timeseries + fits
                ecog_plotPRFtsfits(data2fit, stimulus, results, channels);
                for f = 1:length(f_ind)
                    if length(f_ind) == 1
                        figureName = sprintf('%s_prftimecoursefits', subject);
                    else
                        figureName = sprintf('%s_prftimecoursefits_%d', subject, f);
                    end                
                    set(f, 'Name', figureName);
                    set(f, 'Position', get(0,'screensize'));
                    set(findall(f,'-property','FontSize'),'FontSize',14)
                    saveas(f, fullfile(plotSaveDir, figureName), 'png'); 
                end
                close all;
                
                % PRFs
                coloropt = 0;
                ecog_plotPRFs(results, stimulus, channels, [], [], coloropt)  
                for f = 1:length(f_ind)
                    if length(f_ind) == 1
                        figureName = sprintf('%s_prfs', subject);
                    else
                        figureName = sprintf('%s_prfs_%d', subject, f);
                    end                
                    set(f, 'Name', figureName);
                    set(f, 'Position', get(0,'screensize'));
                    set(findall(f,'-property','FontSize'),'FontSize',14)
                    saveas(f, fullfile(plotSaveDir, figureName), 'png'); 
                end
                close all;
            end
        end
    end
end

end