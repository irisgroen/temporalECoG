function tde_computePRFs(recomputeData, doPlots)

% [results] = tde_computePRFs(recomputeData, recomputeFits, doPlots) 

% <saveDir> path to save parameters and fits; if empty, results are not
%   saved (default)
% <saveName> string to add to the save filename, if results are saved
%   (default empty)

% <recomputeData>
if ~exist('recomputeData','var') || isempty(recomputeData)
    recomputeData = false; % boolean
end

% <doPlots>
if ~exist('doPlots','var') || isempty(doPlots)
    doPlots = false; % boolean
end

%% Set paths
if ~exist('saveDir', 'var'), saveDir = fullfile(analysisRootPath, 'prfs'); end
if ~exist('resultsStr', 'var'), resultsStr = 'prfs'; end
if doPlots
    plotSaveDir = fullfile(analysisRootPath, 'figures', 'prfs');
    if ~exist(plotSaveDir, 'dir'); mkdir(fullfile(plotSaveDir, 'data'));  mkdir(fullfile(plotSaveDir, 'modelfits'));end
end

%% Get the PRF data
recomputeFlag = recomputeData;
subjects      = []; % will default to all subjects in subjectList.tsv
sessions      = []; % will default to all sessions per subject
tasks         = {'prf'};
description   = []; % will default to broadband
epochTime     = [-0.2 0.6];
sampleRate    = []; % will default to 512
saveStr       = 'prfdata';

[fulldata] = tde_getData(recomputeFlag, subjects, sessions, tasks, description, epochTime, sampleRate, saveStr);

% Load the stimulus apertures
stimName = fullfile(tdeRootPath, 'prf_apertures', 'bar_apertures.mat');
load(stimName, 'bar_apertures');
bar_apertures = imresize(bar_apertures, [100 100], 'nearest');

% Define a time window over which to average the broadband timecourse
time_win  = [0.05 0.55];

% Compute PRF timecourses for each subject
nSubjects = length(fulldata);

data2fit = cell(nSubjects,1);
 
for ii = 1:nSubjects
    
    subject  = fulldata{ii}.subject;
    t        = fulldata{ii}.t;
    epochs   = fulldata{ii}.epochs;
    
    % Determine which stimuli to select for the PRF timecourse
    switch subject
        case 'beilen'
            % in this subject, trigger 41 was not sent because it was the
            % same as used for the blank (stimulus coding error).
            stimInx = setdiff(1:224, 41);
        case 'som661'
            % in this subject, more triggers were sent than the actual prf
            % bar positions --> need to interpolate? skip for now
            continue
        otherwise
            stimInx = 1:224;
    end
    
    % Normalize epochs, within run
    [~, nTrials, nChans] = size(epochs);
    nStim = length(stimInx);    
    nRuns = nTrials/nStim; 
    run_indices = []; 
    for jj = 1:nRuns
        run_indices = [run_indices; ones(nStim,1)*jj];
    end
    [epochs] = ecog_normalizeEpochs(epochs, t, [], [], run_indices);
    
    % Compute average broadband response in time window
    trials = squeeze(mean(epochs(t>time_win(1) & t<time_win(2),:,:),1)); 

    % Transpose to have channels in first dimension
    trials = trials';
    
    % Reshape to separate individual runs    
    ts = reshape(trials,[nChans nStim nRuns]);
    
    data2fit{ii}.ts       = ts;
    data2fit{ii}.stim_inx = stimInx;

    % Make plots of the trials and of the PRF timecourses
    if doPlots
        channels = fulldata{ii}.channels;
        events   = fulldata{ii}.events;
        nEpochs = height(events);

        % Plot individual trials
        figureName = sprintf('%s_prftrials', subject);
        figure('Name', figureName); plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
        for el = 1:nChans
            subplot(plotDim1,plotDim2,el); hold on
            for kk = 1:nEpochs
                %ecog_plotSingleTimeCourse(t, epochs(:,kk,el));
                plot(t, epochs(:,kk,el)); 
            end
            plot(t, squeeze(mean(epochs(:,kk,el),2)), 'k', 'LineWidth', 3); axis tight;

            yLim = get(gca, 'YLim');
            line([0 0], yLim,'LineStyle', ':', 'Color', 'k');
            line([t(1) t(end)], [0 0],'LineStyle', ':', 'Color', 'k');
            plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        
            title(plotTitle);
            axis tight
            if el == 1; xlabel('Time (s)'); ylabel('Broadband signal change');end
            set(gcf, 'Position', [150 100 1500 1250]);
            set(gca, 'FontSize', 14);
        end
        saveas(gcf, fullfile(plotSaveDir,'data', figureName), 'png'); close;
        
        % Plot PRF timecourses for each run + average
        figureName = sprintf('%s_prftimecourses', subject);
        figure('Name', figureName); plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
        for el = 1:nChans
            subplot(plotDim1,plotDim2,el); hold on
            plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        
            runNames = [];
            for kk = 1:nRuns
                plot(1:nStim, squeeze(ts(el,:,kk)), 'LineWidth', 1);
                runNames = [runNames {sprintf('run %d', kk)}];
            end
            plot(1:nStim, squeeze(mean(ts(el,:,:),2)), 'k', 'LineWidth', 3); axis tight;
            runNames = [runNames {'average'}];
            title(plotTitle);
            if el == 1; xlabel('PRF stimulus (#)'); ylabel('Broadband signal change'); legend(runNames); end
            set(gcf, 'Position', [150 100 1500 1250]);
            set(gca, 'FontSize', 14);
        end
        saveas(gcf, fullfile(plotSaveDir,'data', figureName), 'png'); close;        
    end   
end
    
% Fit the PRF time courses with analyzePRF
tr             = 1;
opt.hrf        = 1;
opt.maxpolydeg = 0;
opt.xvalmode   = 0; 
opt.display    = 'off';

% Loop over subjects, fit data
for ii = 1:nSubjects

    subject = fulldata{ii}.subject;
	channels = fulldata{ii}.channels;
    
    data = data2fit{ii};
    
    if ~isempty(data)
        
        % Average runs
        data = mean(data.ts,3);
    
        % Define stimulus
        stim_inx = data2fit{ii}.stim_inx;
        stimulus = bar_apertures(:,:,stim_inx);
        %stimulus = {bar_apertures,bar_apertures,bar_apertures,bar_apertures};

        results = analyzePRF(stimulus,data,tr,opt);

        % Save fits to results directory
        if ~isempty(saveDir)

            if ~exist(saveDir, 'dir'); mkdir(saveDir); end

            saveName = sprintf('%s_%s', subject, resultsStr);
            saveName = fullfile(saveDir, saveName);
            fprintf('[%s] Saving results to %s \n', mfilename, saveName);

            if exist(sprintf('%s.mat',saveName),'file')
                warning('[%s] Results file already exists! Writing new file with date-time stamp.',mfilename);
                saveName = sprintf('%s_%s', saveName, datestr(now,30));
                fprintf('[%s] Saving results to %s \n', mfilename, saveName);
            end
            save(saveName, 'channels','data','results','tr','opt','subject');  
        end
    
        % Make plots of the estimated PRFs and PRF fits

        if doPlots
            cfactor = 16.6/100;

            % Timeseries + fits
            figureName = sprintf('%s_prfmodelfits', subject);
            figure('Name', figureName); plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
            
            % Define some variables
            res = [100 100];                    % row x column resolution of the stimuli
            resmx = 100;                        % maximum resolution (along any dimension)
            hrf = results.options.hrf;          % HRF that was used in the model
            degs = results.options.maxpolydeg;  % vector of maximum polynomial degrees used in the model

            % Pre-compute cache for faster execution
            [d,xx,yy] = makegaussian2d(resmx,2,2,2,2);

            % Prepare the stimuli for use in the model
            stimulusPP = {};
            for cc=1:length(stimulus)
              stimulusPP{cc} = squish(stimulus{cc},2)';  % this flattens the image so that the dimensionality is now frames x pixels
              stimulusPP{cc} = [stimulusPP{cc} cc*ones(size(stimulusPP{cc},1),1)];  % this adds a dummy column to indicate run breaks
            end

            % Define the model function.  This function takes parameters and stimuli as input and
            % returns a predicted time-series as output.  Specifically, the variable <pp> is a vector
            % of parameter values (1 x 5) and the variable <dd> is a matrix with the stimuli (frames x pixels).
            % Although it looks complex, what the function does is pretty straightforward: construct a
            % 2D Gaussian, crop it to <res>, compute the dot-product between the stimuli and the
            % Gaussian, raise the result to an exponent, and then convolve the result with the HRF,
            % taking care to not bleed over run boundaries.
            modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));

            % Construct projection matrices that fit and remove the polynomials.
            % Note that a separate projection matrix is constructed for each run.
            polymatrix = {};
            for cc=1:length(degs)
              polymatrix{cc} = projectionmatrix(constructpolynomialmatrix(size(data{cc},2),0:degs(cc)));
            end

            for el = 1:nChans
                subplot(plotDim1,plotDim2,el); hold on
                plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        

                vx = el;
                % For each run, collect the data and the model fit.  We project out polynomials
                % from both the data and the model fit.  This deals with the problem of
                % slow trends in the data.
                datats = {};
                modelts = {};
                for cc=1:length(data)
                  datats{cc} =  polymatrix{cc}*data{cc}(vx,:)';
                  modelts{cc} = polymatrix{cc}*modelfun(results.params(1,:,vx),stimulusPP{cc});
                end

                set(gcf,'Units','points');
                plot(cat(1,datats{:}),'k-', 'LineWidth', 2);
                plot(cat(1,modelts{:}),'r-','LineWidth', 2);
                xlabel('PRF stimulus','FontSize', 28);
                if el == 1, legend('Data', 'Model prediction'); end
                set(gca, 'FontSize', 18)
                title(plotTitle);
                axis tight
            end
            saveas(gcf, fullfile(plotSaveDir,'modelfits', figureName), 'png'); close;

            % PRFs
            coloropt = 0;
            figureName = sprintf('%s_prfs', subject);
            figure('Name', figureName); plotDim1 = round(sqrt(nChans)); plotDim2 = ceil((nChans)/plotDim1);
            for el = 1:nChans
                subplot(plotDim1,plotDim2,el); hold on
                plotTitle = sprintf('%s %s %s ', channels.name{el}, channels.bensonarea{el}, channels.wangarea{el});        

                sd = results.rfsize(el);
                [xx, yy] = meshgrid(linspace(-1,1,100));
                [th, r] = cart2pol(xx, yy);
                p = results.params(1,:,el);
                im = makegaussian2d(250,p(1)+75,p(2)+75,p(3)/sqrt(p(5)),p(3)/sqrt(p(5))); 

                % plot pRF
                imagesc(im);
                if coloropt == 1
                    colormap(parula)
                else
                    colormap(1-gray)
                end

                % draw stimulus
                h1 = k_drawellipse(125,125,0,50,50); % circle indicating stimulus extent
                set(h1,'Color',[0 0 0],'LineWidth',1, 'LineStyle', ':');
                h2 = straightline(125,'h','k:');     % line indicating horizontal meridian
                h3 = straightline(125,'v','k:');     % line indicating vertical meridian
                if coloropt == 1
                    set(h1,'Color',[1 1 1]);
                    set(h2,'Color',[1 1 1]);
                    set(h3,'Color',[1 1 1]);
                end
                % plot pRF center and sd  
                h1 = k_drawellipse(p(2)+75,p(1)+75,0,sd,sd);      % 
                h2 = k_drawellipse(p(2)+75,p(1)+75,0,2*sd,2*sd);  % 
                set(h1,'Color', [0 0 0],'LineWidth',2,'LineStyle', '-');
                set(h2,'Color', [0 0 0],'LineWidth',2,'LineStyle', '-');
                h3 = scatter(p(2)+75,p(1)+75,'wo','filled');
                if coloropt == 1
                    set(h3,'CData',[1 0 0]);
                end

                axis square;
                set(gca, 'XTick', [1 25:25:250], 'XTickLabel', round(([1 25:25:250]*0.16)-20));
                set(gca, 'YTick', [1 25:25:250], 'YTickLabel', round(([1 25:25:250]*0.16)-20));
                xlabel('X-position (deg)');
                ylabel('Y-position (deg)');
                title(plotTitle);
                set(gca, 'XLim', [25 225],'YLim',[25 225])
            end
            saveas(gcf, fullfile(plotSaveDir,'modelfits', figureName), 'png'); close;
            
            % Add some fun summary plots compaing e.g. 
        end
    end
end