function [chan_idx, channels, group_prob] = groupElecsByVisualArea(channels, groupingMethod, areaNames, nResamples)
% Groups electrodes in a channel table based on atlas assignment
%
% [chan_idx, channels, group_prob] = groupElecsByVisualArea(channels,
%                                       [groupingMethod], [areaNames])
%
% Input
%     channels:         BIDSformatted channels table (matlab table)
%     groupingMethod:   electrode grouping method (string). options
%                           - 'fixedassignment' (default) requires
%                               wang_15_mplbl and benson14_varea column 
%                               in channel table
%                           - 'probabilisticresample'; requires
%                               wang_15_fplbl column in channel table
%     areaNames:        list of areas to group (cell array of string)
%                               default: {'V1', 'V2', 'V3', 'V3a','V3b',...
%                                         'LO1','LO2','TO','IPS'};
%
% Output
%     chan_idx:         for 'fixedassignment':
%                               channels * areas mask (logicals) indicating
%                               which channels to include in which areas
%                       for 'probabilistic resample':
%                               channels * areas * nResamples mask indicating
%                               which channels to include in which areas
%                               for nResamples resamples of the channels
%     channels:         updated channels table with a single row per area
%                       instead of individual channels, and new columns
%                       'subject_name' (indicating which subjects were
%                       included in which area) and 'number_of_elecs'
%                       (indicating how many electrodes were included).
%                       NOTE: for 'probabilisticresample', the subject_name
%                       and nelecs columns contain ANY subject and the 
%                       MAXIMUM number of included elecs in each area
%     group_prob:       probability of belonging to each of the areaNames
%                       according to the Wang 2015 full probability atlas. Note
%                       that the probalities here are normalized such that
%                       the probability of belonging to no area is removed,
%                       and that some areas are grouped (e.g. V1d+V1v =
%                       V1); see matchAreaNameToAtlas.m
%
% IG 2020

if ~exist('groupingMethod', 'var') || isempty(groupingMethod)
    groupingMethod = 'fixedassignment';
end

if ~exist('areaNames', 'var') || isempty(areaNames)
    areaNames = {'V1','V2','V3','V3a','V3b','LO1','LO2','TO','IPS'};
end

if ~exist('nResamples', 'var') || isempty(nResamples)
    nResamples = 1000;
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
        
        % For each atlas, determine which ROIs match the requested areaNames
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
        group_prob = [];
        
    case 'probabilisticresample' 
        
        atlasName = 'wang15_fplbl';
        
        % Check that the requested atlas column is present
        assert(any(contains(channels.Properties.VariableNames, atlasName)));
        
        % Get the wang probability values
        idx = contains(channels.Properties.VariableNames, atlasName);

        % Normalize probabilities by the sum (i.e. exclude 'none')
        probvals = channels(:,idx).Variables;
        probvals_norm = probvals./sum(probvals,2);
        probvals_norm(isnan(probvals_norm)) = 0;
        
        % Group wang ROIs into bigger groups (e.g. V1v + V1d -> V1)
        atlasLabels = getAtlasLabels(atlasName);        
        probvals_norm_grouped = nan(height(channels), length(areaNames));     
        for ii = 1:nAreas
            area_idx = matchAreaNameToAtlas(areaNames{ii}, atlasLabels);
            probvals_norm_grouped(:,ii) = sum(probvals_norm(:,area_idx),2);
        end
        
        % Group wang ROIs into bigger groups (e.g. V1v + V1d -> V1)
        group_prob = probvals_norm_grouped;

        % Generate nResamples channel indices for each area
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
        % channels table below will contain ANY subject included in any
        % area, and the MAXIMUM number of included elecs in each area
       
end

% Track which subjects were included in which area
subjects = cell(nAreas,1); 
nelecs = nan(nAreas,1); 
for ii = 1:nAreas
    temp = unique(channels.subject_name(chan_idx_mask(:,ii)));
    subjects{ii} = [temp{:}];
    nelecs(ii) = length(find(chan_idx_mask(:,ii)));
end

% Create a new channels table:
name               = areaNames';
type               = repmat({'n/a'}, [nAreas 1]);
units              = repmat(channels.units(1), [nAreas 1]);
sampling_frequency = repmat(channels.sampling_frequency(1), [nAreas 1]);
subject_name       = subjects;    
number_of_elecs    = nelecs;

channels = table(name, type, units, sampling_frequency, subject_name, number_of_elecs);
index = [1:height(channels)]';
channels = addvars(channels, index, 'Before', 'name'); 

end

        