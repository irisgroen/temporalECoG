function [epochs, channels, epochs_se] = tde_prepareData(epochs, channels, opts)
% Now that we have the selected trials and electrodes from all subjects, 
% run a few checks, and make decisions on how to further process the data.
%
% For example:
% Sort electrodes on area (optional)
% Normalizes by max (optional)
% Average across visual areas (optional)

if ~isfield(opts,'normalize_data') || isempty(opts.normalize_data)
    opts.normalize_data      = false;  % boolean
end
if ~isfield(opts,'average_elecs') || isempty(opts.average_elecs)
    opts.average_elecs       = false; % boolean
end
if ~isfield(opts,'sort_channels') || isempty(opts.sort_channels)
    opts.sort_channels       = true;  % boolean
end
if ~isfield(opts,'area_names') || isempty(opts.area_names)
    area_names = [];
end

% Sort electrodes on visual area (rather than subjectID)?
if opts.sort_channels
    % sort on benson area
    [~,I] = sortVisualAreaNames(channels.benson14_varea);
    channels = channels(I,:);
    epochs = epochs(:,:,I);
end

% Average elecs within area?
if opts.average_elecs
    fprintf('[%s] Averaging electrodes...\n', mfilename);
    [epochs, channels, epochs_se] = average_elecs(epochs, channels, area_names);
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
        tmp = data(:, :, ii);
        maxRsp(ii) = max(tmp(:));
        normdata(:, :, ii) = data(:, :, ii)./maxRsp(ii);
    end
    data = normdata;
end

% Data averaging
function [data, channels, se] = average_elecs(data, channels, area_names)
    
    [~, channels, group_prob] = groupElecsByVisualArea(channels, 'probabilisticresample', area_names);
    fun = @mean;
    [data, se] = averageWithinArea(data, group_prob, fun);   
end
