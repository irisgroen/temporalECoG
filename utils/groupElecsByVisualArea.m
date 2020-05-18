function [chan_idx, channels, areaNames] = groupElecsByVisualArea(channels, groupingMethod, areaNames)
% Groups electrodes in a channel table based on atlas values
%
% [chan_idx, channels, areaNames] = groupElecsByVisualArea(channels,
% [groupingMethod], [areaNames])
%
% Input
%     channels:         BIDSformatted channels table (matlab table)
%     groupingMethod:   electrode grouping method (string). options
%                           - 'fixedassignment'; requires wang_15_mplbl and
%                               benson14_varea column in channel table
%                           - 'probabilisticresample'; requires
%                               wang_15_fplbl column in channel table
%                               default: 'fixedassignment'
%     areaNames:        list of areas to group (cell array of string)
%                               default: {'V1', 'V2', 'V3', 'V3a','V3b',...
%                                         'LO1','LO2','TO1','IPS'};
%
% Input
%     chan_idx:         for 'fixedassignment':
%                               channels * areas mask (logicals) indicating
%                               which channels to include in which areas
%                       for 'probabilistic resample':
%                               channels * areas * 100 mask indicating
%                               which channels to include in which areas
%                               for 100 resamples of the channels
%     channels:         updated channels table with a single row per area
%                       instead of individual channels, and new columns
%                       'subject_name' (indicating which subjects were
%                       included in which area) and 'number_of_elecs'
%                       (indicating how many electrodes were included
%                       NOTE: for 'probabilisticresample', the subject_name
%                       and nelecs columns contain ANY subject included in
%                       any area, and the MAXIMUM number of included elecs
%                       in each area
%     areaNames:        list of areas used for grouping 
%
% IG 2020

if ~exist('groupingMethod', 'var') || isempty(groupingMethod)
    groupingMethod = 'fixedassignment';
end

if ~exist('areaNames', 'var') || isempty(areaNames)
    areaNames = {'V1', 'V2', 'V3', 'V3a', 'V3b','LO1','LO2','TO1','IPS'};
end

chan_idx = [];
nAreas = length(areaNames);

switch groupingMethod
    
    case 'fixedassignment' % for each area, find matched electrodes according to either wang15_mplpl or benson14_varea
        
        atlasNames = {'wang15_mplbl', 'benson14_varea'};
        nAtlases = length(atlasNames);
        
        % Check that the requested columns are present
        for ii = 1:nAtlases
            assert(any(contains(channels.Properties.VariableNames, atlasNames{ii})))
        end
        
        % For each atlas, determine with ROIs match the requested areaNames
        for ii = 1:nAtlases
            atlasLabels{ii} = getAtlasLabels(atlasNames{ii});
            for jj = 1:nAreas
                area_idx{ii}{jj} = matchAreaNameToAtlas(areaNames{jj}, atlasLabels{ii});
            end
        end
        
        % Match the requested areaNames to the channel area labels
        for jj = 1:nAreas
            for ii = 1:nAtlases
                matched_idx(:,ii) = matchAreaNameToAtlas(atlasLabels{ii}(area_idx{ii}{jj}),channels.(atlasNames{ii}));
            end
            chan_idx(:,jj) = any(matched_idx,2);
        end
        chan_idx = logical(chan_idx);
        chan_idx_mask = chan_idx;
        
    case 'probabilisticresample' 
        
        atlasName = 'wang15_fplbl';
        
        % Check that the requested atlast column is present
        assert(any(contains(channels.Properties.VariableNames, atlasName)));
        
        % Get the wang probability values
        idx = contains(channels.Properties.VariableNames, atlasName);

        % Normalize probabilities by the sum (i.e. exclude 'none')
        probvals = channels(:,idx).Variables;
        probvals_norm = probvals./sum(probvals,2);
        
        % Group wang ROIs into bigger groups (e.g. V1v + V1d -> V1)
        atlasLabels = getAtlasLabels(atlasName);        
        probvals_norm_grouped = nan(height(channels), length(areaNames));     
        for ii = 1:nAreas
            area_idx = matchAreaNameToAtlas(areaNames{ii}, atlasLabels);
            probvals_norm_grouped(:,ii) = sum(probvals_norm(:,area_idx),2);
        end
          
        % Generate nResamples channel indices for each area
        nResamples = 100;
        probn = round(probvals_norm_grouped*nResamples);
        
        % Generate a set of resampled assignments
        nChans = size(probn,1);
        nAreas = size(probn,2);
        chan_idx = zeros([nChans nAreas nResamples]);
        
        for ii = 1:nChans
            for jj = 1:nAreas
                inx = randperm(nResamples, probn(ii,jj));
                chan_idx(ii,jj,inx) = 1;
            end
        end
        chan_idx = logical(chan_idx);
        
        % NOTE: I tried to do the resampled assignment without a forloop
        % but couldn't get it to work
%         k = nResamples;
%         x = probn; % also tried prob(probn>0)
%         B = arrayfun(@(z) randsample(k, z), x, 'UniformOutput', 'false');
        
        chan_idx_mask = probn>0;        
        % NOTE: this means that the subjects and nelecs columns in the
        % channels table below will contain ANY subject include in any
        % area, and the MAXIMUM number of included elecs in each area
       
end

% Track which subjects were included in which area
subjects = cell(nAreas,1); 
nelecs = cell(nAreas,1); 
for ii = 1:nAreas
    temp = unique(channels.subject_name(chan_idx_mask(:,ii)));
    subjects{ii} = [temp{:}];
    nelecs{ii} = length(find(chan_idx_mask(:,ii)));
end

% Create a new channels table:
name               = areaNames';
type               = repmat({'n/a'}, [nAreas 1]);
units              = repmat(channels.units(1), [nAreas 1]);
sampling_frequency = repmat(channels.sampling_frequency(1), [nAreas 1]);
subject_name       = subjects;    
number_of_elecs    = nelecs;

channels = table(name, type, units, sampling_frequency, subject_name, number_of_elecs);

end

%% SUBROUTINES

function area_idx = matchAreaNameToAtlas(areaName, atlasLabels)
    % find atlasLabels for that contain a given areaName (e.g. V1v and V1d
    % for V1 for wang atlases), preventing matches for V3a/b with V3 and
    % collapsing all IPS maps
    if strcmpi(areaName, 'V3')
        area_idx = contains(atlasLabels, 'V3') & ~contains(atlasLabels, {'V3a', 'V3b'});
    elseif strcmpi(areaName, 'IPS')
        area_idx = contains(atlasLabels, {'IPS0', 'IPS1', 'IPS2', 'IPS3', 'IPS4', 'IPS5'});
    else
        area_idx = contains(atlasLabels, areaName); 
    end
end
        