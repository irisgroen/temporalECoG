function [m, se, indiv_points] = averageWithinArea(data, chan_idx, fun)
    
% data: e.g. params * channels (last dim should be channels)
% chan_idx: as outputted by groupElecsByVisualArea.m

if ~exist('fun','var') || isempty(fun)
    fun = @median;
end

if ndims(data) == 3
    multiDimData = true;
    dataSz = size(data);
    data = reshape(data, [dataSz(1) * dataSz(2) dataSz(3)]);
else
    multiDimData = false;
end

[~, nAreas, nResamples] = size(chan_idx);

m = nan(size(data,1), nAreas); 
se = nan(size(data,1), nAreas, 2);
indiv_points = cell(nAreas,1);

if nResamples == 1
    
    % bootstrap across electrodes
    numboot = 1000;
    for ii = 1:nAreas
        elec_idx = chan_idx(:,ii);
        %[m(:,ii),se(:,ii,:)] = calcmdsepct(data(:,elec_index),2);
        d = data(:,elec_idx);
        m(:,ii) = median(d,2);
        func = @(x) median(x,2,'omitnan');
        mdata = bootstrp(numboot,func,d)';
        se(:,ii,:) = squeeze(prctile(mdata, [15.87 84.13],2));   
        indiv_points{ii} = data(:,elec_idx);
    end
        
else
    
    % resample use based on wang probabilities (new method)
        
    for ii = 1:nAreas
        fprintf('[%s] computing probabilistic average for area %d \n', mfilename, ii);
        % get index of channels for this area
        inx = squeeze(chan_idx(:,ii,:)); 
        % generate nresample copies of the data
        data_resampled = repmat(data, [ones(1, ndims(data)) nResamples]);
        % set non-included resamples to nan
        data_resampled(:,~inx) = nan;
        % take mean/median across the electrodes nresample times        
        mdata = squeeze(fun(data_resampled,2,'omitnan'));
        % take median and 68% confidence interval of resampled
        % distribution of electrode means
        m(:,ii) = median(mdata,2, 'omitnan');
        se(:,ii,:) = prctile(mdata,[15.87 84.13],2);         
        % include any included electrode in the individual points
        elec_idx = any(inx,2);
        indiv_points{ii} = data(:,elec_idx);
    end    
end

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
