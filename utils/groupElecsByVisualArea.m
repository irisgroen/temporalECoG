function [INX, channels] = groupElecsByVisualArea(channels, areanames)

% ADD FUNCTION DESCRIPTION

if ~exist('areanames', 'var') || isempty(areanames)
    areanames = {'V1', 'V2', 'V3', 'V3a', 'V3b','LO1','LO2','TO1','IPS'};
end

INX = [];
for ii = 1:length(areanames)
    if strcmpi(areanames{ii}, 'V3')
        INX{ii} = contains(channels.wangarea, 'V3') & ~contains(channels.bensonarea, {'V3a', 'V3b'}) | contains(channels.bensonarea, 'V3') & ~contains(channels.bensonarea, {'V3a', 'V3b'});
    elseif strcmpi(areanames{ii}, 'IPS')
        INX{ii} = contains(channels.wangarea, {'IPS0', 'IPS1', 'IPS2', 'IPS3', 'IPS4', 'IPS5'});
    else
        INX{ii} = contains(channels.wangarea, areanames{ii}) | contains(channels.bensonarea, areanames{ii});
    end
end
    
subjects = cell(length(INX),1); 
nelecs = cell(length(INX),1); 

for ii = 1:length(INX)
    temp = unique(channels.subject_name(INX{ii}));
    subjects{ii} = [temp{:}];
    nelecs{ii} = length(find(INX{ii}));
end

% Create a new channels table:
name               = areanames';
type               = repmat({'n/a'}, [length(name) 1]);
units              = repmat(channels.units(1), [length(name) 1]);
sampling_frequency = repmat(channels.sampling_frequency(1), [length(name) 1]);
subject_name       = subjects;    
number_of_elecs    = nelecs;

channels = table(name, type, units, sampling_frequency, subject_name, number_of_elecs);

end