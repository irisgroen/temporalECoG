function [m_boot, se] = averageWithinArea(data, group_prob, fun, numboot)
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
%     numboot:          (optional) Number of bootstraps (default: 10000).
%
% Output
%     m:                estimated distribution metric (default median).
%     se:               68% confidence interval.
%
% IG 2020

if ~exist('fun','var') || isempty(fun)
    fun = @median;
end

if ~exist('numboot','var') || isempty(numboot)
    numboot = 10000;
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

% each boot, sample 112 at once
for ii = 1:numboot
    elec_idx = randsample(1:nElecs,nElecs,1);
    sampled_data = data(:,elec_idx);
    % assign to area probabilistically
    elec_area = assign_to_area(group_prob(elec_idx,:)); 
    % do summary per area - store result
    for jj = 1:nAreas
        md(:,jj,ii) = fun(sampled_data(:,elec_area == jj),2);
    end
end

m_boot = fun(md,3, 'omitnan');
se = squeeze(prctile(md,[15.87 84.13],3));

% reshape back to original shape
if multiDimData
    m_boot = reshape(m_boot, [dataSz(1) dataSz(2) nAreas]); 
    se = reshape(se, [dataSz(1) dataSz(2) 2 nAreas]);
end

end

function elec_area = assign_to_area(group_prob)
    % sample a random number between 0 and 1 
    p = rand(size(group_prob,1),1);% * ones(1,size(gp,2));
    % compute cumulative sum
    a = cumsum(group_prob,2) > p;
    sa = sum(a,2) < 1;
    [~,elec_area] = max(a,[],2);
    elec_area(sa) = nan;
end
