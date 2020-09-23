function [m, se] = averageWithinArea(data, group_prob, fun)
% Averages data across electrodes belonging to individual areas. The
% assignment to areas is based on the output and method specified in groupElecsByVisualArea.m
%
% [m, se] = averageWithinArea(data, group_prob, [fun])
%
% Input
%     data:             The to be averaged data, e.g. a list of fitted
%                       parameters for each channel. Can be multi-
%                       dimensional (e.g. a time-course for each channel),
%                       but the last dimension should be channels.
%     group_prob:       Matrix of size elec x areas containing electrode
%                       probabilities of belonging to a set of areas as
%                       outputted by groupElecsByVisualArea.m
%     fun:              (optional) Averaging function (default: @median).
%
% Output
%     m:                estimated area average (default median).
%     se:               68% confidence interval.
%
% IG 2020

if ~exist('fun','var') || isempty(fun)
    fun = @median;
end

% if there are more than one dimension per channel, vectorize data
if ndims(data) == 3
    multiDimData = true;
    dataSz = size(data);
    data = reshape(data, [dataSz(1) * dataSz(2) dataSz(3)]);
else
    multiDimData = false;
end

% pre-allocate
[~, nAreas] = size(group_prob);
[~, nElecs] = size(data);
m = nan(size(data,1), nAreas); 
se = nan(size(data,1), nAreas, 2);

% generate averages
numsamp = 100;
numboot = 1000;

for ii = 1:nAreas
    % sample with replacement from the total number electrodes, with a
    % weight based on their probability of belonging to this area
    elec_idx = randsample(1:nElecs,numsamp,1,group_prob(:,ii));
    data_sample = data(:,elec_idx)'; % put channel dimension first
    % take median across the probabilistically sampled channels for this area
    m(:,ii) = fun(data_sample);
    % generate a distribution of medians through bootstrapping
    mdata = bootstrp(numboot,fun,data_sample)';
	% take 68% confidence interval of resampled distribution of medians
	%se(:,ii,:) = squeeze(prctile(mdata,[15.87 84.13],2));
    % take 95% confidence interval of resampled distribution of medians
    se(:,ii,:) = squeeze(prctile(mdata, [0.05 0.95],2));
end

% if nResamples == 1
%     
%     % bootstrap across electrodes (fixed assignment)
%     
%     for ii = 1:nAreas
%         fprintf('[%s] computing bootstrapped average based on fixed assignment for area %d \n', mfilename, ii);
%         elec_idx = chan_idx(:,ii);  
%         d = data(:,elec_idx); 
%         d = d'; % make sure the channel dimension is first
%         % take true median across the channels for this area
%         m(:,ii) = fun(d)'; 
%         % generate a distribution of medians through bootstrapping
%         mdata = bootstrp(numboot,fun,d)';
%         % take 68% confidence interval of resampled
%         % distribution of electrode means
%         se(:,ii,:) = squeeze(prctile(mdata, [15.87 84.13],2));
%         % also return the individual channel values
%         indiv_points{ii} = data(:,elec_idx);
%         
%          % debug:
%         % I compared this with knkutils function:
%         %   [m(:,ii),se(:,ii,:)] = calcmdsepct(data(:,elec_index),2);
%         % which gives the same result but uses parpool which is slow to
%         % start up.
%     end
%         
% else
%     
%     % resample based on wang probabilities (probabilistic resample)
%         
%     for ii = 1:nAreas
%         
%         fprintf('[%s] computing bootstrapped average based on probabilistic assignment for area %d \n', mfilename, ii);
%         m_resamp = [];
%         mdata_resamp = [];
%         for jj = 1:nResamples
%             % get index of channels for this area and resample 
%             elec_idx = squeeze(chan_idx(:,ii,jj));
%             % select the data for these electrodes
%             d = data(:,elec_idx);
%             d = d'; % make sure the channel dimension is first
%             % take true median across the channels for this resample
%             m_resamp(:,jj) = fun(d); 
%             % generate a distribution of medians through bootstrapping
%             mdata_resamp(:,:,jj) = bootstrp(numboot,fun,d)';         
%         end
%         % take true median across the resamples for this area
%         m(:,ii) = fun(m_resamp');
%         % collapse resamples and bootstraps into one distribution
%         mdata = mdata_resamp(:,:);
%         % take 68% confidence interval of resampled distribution            
%         se(:,ii,:) = squeeze(prctile(mdata, [15.87 84.13],2));
%         % also return the individual channel values
%         indiv_points{ii} = data(:,elec_idx);
%       
% %         fun = @(x) median(x,2,'omitnan');
% % 
% %         %get index of channels for this area 
% %         inx = squeeze(chan_idx(:,ii,:));  
% %         % generate nresample copies of the data
% %         data_resampled = repmat(data, [ones(1, ndims(data)) nResamples]);
% %         % set non-included resamples to nan
% %         data_resampled(:,~inx) = nan;
% %         % take median across the electrodes nresample times        
% %         mdata = squeeze(fun(data_resampled));
% %        
% %         if size(mdata,2) == 1, mdata = mdata'; end     
% %         % take median and 68% confidence interval of resampled
% %         % distribution of electrode means
% %         m(:,ii) = fun(mdata);
% %         se(:,ii,:) = prctile(mdata,[15.87 84.13],2);         
% %         %include any included electrode in the individual points
% %         elec_idx = any(inx,2);
% %         
% %         indiv_points{ii} = data(:,elec_idx);
%     end    
% end

% reshape back to original shape
if multiDimData
    m = reshape(m, [dataSz(1) dataSz(2) nAreas]); 
    se = reshape(se, [dataSz(1) dataSz(2) 2 nAreas]);
end

% debug

% figure;hold on
% errorbar(1:nAreas, m, m-se(:,:,1), se(:,:,2)-m, '.b', 'MarkerSize', 30, 'LineWidth', 4, 'LineStyle', 'none', 'CapSize', 0)
% errorbar(1:nAreas, m2, m2-se2(:,:,1), se2(:,:,2)-m2, '.c', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
% set(gca, 'XLim', [0 nAreas+1], 'XTick', 1:nAreas, 'XTickLabel', areaNames');
% set(gca, 'YLim', [0 1]);
% legend({'probresample', 'oldmethod'})

end
