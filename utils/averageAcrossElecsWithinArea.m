function [m, se, indiv_points] = averageAcrossElecsWithinArea(data, chan_idx)
    
% data: e.g. params * channels (last dim should be channels)
% chan_idx: as outputted by groupElecsByVisualArea.m


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
        m(:,ii) = median(d);
        func = @(x) median(x,'omitnan');
        mdata = bootstrp(numboot,func,d);
        se(:,ii,:) = prctile(mdata, [15.87 84.13]);   
        indiv_points{ii} = data(:,elec_idx);

    end
        
else
    
    % resample use based on wang probabilities (new method)
    
    % generate nresample copies of the data 
    data_resampled = repmat(data', [1 nResamples]);
       
    for ii = 1:nAreas
        inx = squeeze(chan_idx(:,ii,:)); 
        % set non-included resamples to nan
        data_resampled(~inx) = nan;
        % take mean across the electrodes nresample times
        mdata = median(data_resampled,'omitnan');
        % take median and 68% confidence interval of resampled
        % distribution of electrode means
        m(:,ii) = median(mdata, 'omitnan');
        se(:,ii,:) = prctile(mdata, [15.87 84.13]);  
        % include any included electrode in the individual points
        elec_idx = any(inx,2);
        indiv_points{ii} = data(:,elec_idx);
    end    
end

% debug

% figure;hold on
% errorbar(1:nAreas, m, m-se(:,:,1), se(:,:,2)-m, '.b', 'MarkerSize', 30, 'LineWidth', 4, 'LineStyle', 'none', 'CapSize', 0)
% errorbar(1:nAreas, m2, m2-se2(:,:,1), se2(:,:,2)-m2, '.c', 'MarkerSize', 30, 'LineWidth', 2, 'LineStyle', 'none', 'CapSize', 0)
% set(gca, 'XLim', [0 nAreas+1], 'XTick', 1:nAreas, 'XTickLabel', areaNames');
% set(gca, 'YLim', [0 1]);
% legend({'probresample', 'oldmethod'})

end
