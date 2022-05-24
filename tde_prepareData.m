function [epochs, channels, epochs_se] = tde_prepareData(epochs, channels, opts)

% Now that we have the selected trials and electrodes from all subjects, 
% run a few checks, and make decisions on how to further process the data.
%
% For example:
% - Sort electrodes on area (opts.sort_channels)
% - Normalize by max (opts.normalize_data)
% - Average across visual areas (opts.average_elecs)
% Note: Adds additional index column to list sorted elecs
% 
% 2022 Iris Groen

if ~isfield(opts,'sort_channels') || isempty(opts.sort_channels)
    opts.sort_channels       = true;  % boolean
end
if ~isfield(opts,'normalize_data') || isempty(opts.normalize_data)
    opts.normalize_data      = false;  % boolean
end
if ~isfield(opts,'average_elecs') || isempty(opts.average_elecs)
    opts.average_elecs       = false; % boolean
end

% Sort electrodes on visual area?
if opts.sort_channels
    % sort on benson area
    [~,I] = sortVisualAreaNames(channels.benson14_varea);
    channels = channels(I,:);
    epochs = epochs(:,:,I);
end

% Average elecs within area?
if opts.average_elecs
    fprintf('[%s] Averaging electrodes...\n', mfilename);
    [epochs, channels, epochs_se] = average_elecs(epochs, channels, opts.areanames);
else
    epochs_se = [];
end

% Scale each electrode to its max?
if opts.normalize_data       
    [epochs] = normalize_data(epochs);
end

% Add an index column to channels
index = [1:height(channels)]';
channels = addvars(channels, index, 'Before', 'name'); 

fprintf('[%s] Done! \n',mfilename);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% SUBROUTINES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% Data normalization
function [data] = normalize_data(data)
    normdata = data;
    % normalize each channel separately
    for ii = 1:size(data,3)
        tmp = data(:,:,ii);
        %normdata(:,:,ii) = data(:,:,ii)./max(tmp(:));
        normdata(:,:,ii) = data(:,:,ii)./norm(tmp(:),2);
    end
    data = normdata;
end

% Data averaging
function [data, channels, se] = average_elecs(data, channels, area_names)
    
    [~, channels, group_prob] = groupElecsByVisualArea(channels, 'probabilisticresample', area_names);
    fun = @mean;
    [data, se] = averageWithinArea(data, group_prob, fun);   
end
