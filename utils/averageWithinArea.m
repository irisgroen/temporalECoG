function [m, se, m_boot, n_elecs_selected] = averageWithinArea(data, group_prob, fun, numboot)
% Averages data across electrodes belonging to individual areas. The
% assignment to areas is based on probability of belonging to each area.
%
% [m, se] = averageWithinArea(data, group_prob, [fun], [numboot])
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
%     numboot:          (optional) Number of bootstraps (default: 1000).
%
% Output
%     m:                estimated summary metric (default median).
%     se:               68% confidence interval.
%
% IG 2020

if ~exist('fun','var') || isempty(fun)
    fun = @median;
end

if ~exist('numboot','var') || isempty(numboot)
    numboot = 1000;
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
m_boot = nan(size(data,1), nAreas, numboot); 
n_elecs_selected = nan(nAreas, numboot);

if numboot > 1
    fprintf('[%s] Computing %s using %d bootstraps...\n', mfilename, func2str(fun), numboot);

    % each boot, sample all channels at once
    for ii = 1:numboot
        elec_idx = randsample(1:nElecs,nElecs,1);
        sampled_data = data(:,elec_idx);
        % assign to area probabilistically
        elec_area = assignElecToAreaProb(group_prob(elec_idx,:)); 
        % do summary per area - store result
        for jj = 1:nAreas
            m_boot(:,jj,ii) = fun(sampled_data(:,elec_area == jj),2);
            % track how many electrodes were sampled on each bootstrap
            n_elecs_selected(jj,ii) = length(find(elec_area == jj));
        end
    end
    
    % take mean of bootstrapped distribution and compute confidence intervals
    m = mean(m_boot,3, 'omitnan');
    se = squeeze(prctile(m_boot,[15.87 84.13],3));
    
    % report average number of electrodes sampled per bootstrap
    for jj = 1:nAreas
        mn_elecs_selected = median(n_elecs_selected(jj,:),2);
        fprintf('[%s] Median no of included elecs in area %d = %0.1f \n',mfilename, jj, mn_elecs_selected);
    end

else
    fprintf('[%s] Computing %s using all electrodes assigned once\n', mfilename, func2str(fun));
    % if not bootstrapping, just use all electrodes and assign once
    m = nan(size(data,1), nAreas);
    while any(isnan(m(:)))
        elec_area = assignElecToAreaProb(group_prob); 
        for jj = 1:nAreas
            m(:,jj) = fun(data(:,elec_area == jj),2, 'omitnan');
            n_elecs_selected(jj) = length(find(elec_area == jj));
            se = [];
        end
    end
end

% reshape back to original shape
if multiDimData
    m = reshape(m, [dataSz(1) dataSz(2) nAreas]); 
    if ~isempty(se)
        se = reshape(se, [dataSz(1) dataSz(2) 2 nAreas]);
    end
end

end
