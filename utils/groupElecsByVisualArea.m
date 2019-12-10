function [INX, channels] = groupElecsByVisualArea(channels)
    

INX = [];
INX{1} = contains(channels.wangarea, 'V1') | contains(channels.bensonarea, 'V1');
INX{2} = contains(channels.wangarea, 'V2') | contains(channels.bensonarea, 'V2');
INX{3} = contains(channels.wangarea, 'V3') & ~contains(channels.bensonarea, {'V3a', 'V3b'}) | contains(channels.bensonarea, 'V3') & ~contains(channels.bensonarea, {'V3a', 'V3b'});
INX{4} = contains(channels.wangarea, {'V3a'}) | contains(channels.bensonarea, {'V3a'}) ;
INX{5} = contains(channels.wangarea, {'V3b'}) | contains(channels.bensonarea, {'V3b'}) ;
INX{6} = contains(channels.wangarea, {'hV4'}) | contains(channels.bensonarea, {'hV4'});
INX{7} = contains(channels.wangarea, {'LO1'}) | contains(channels.bensonarea, {'LO1'});
INX{8} = contains(channels.wangarea, {'LO2'}) | contains(channels.bensonarea, {'LO2'});
INX{9} = contains(channels.wangarea, {'TO1'}) | contains(channels.bensonarea, {'TO1'});
%INX{10} = contains(channels.wangarea,{'TO2'}) | contains(channels.bensonarea, {'TO2'});
INX{10} = contains(channels.wangarea, {'IPS0', 'IPS1', 'IPS2', 'IPS3', 'IPS4', 'IPS5'});
%INX{12} = contains(channels.wangarea, {'VO1','VO2', 'PHC1', 'PHC2'}) | contains(channels.bensonarea, {'VO1', 'VO2'});

subjects = cell(length(INX),1); 
nelecs = cell(length(INX),1); 

for ii = 1:length(INX)
    temp = unique(channels.subject_name(INX{ii}));
    subjects{ii} = [temp{:}];
    nelecs{ii} = length(find(INX{ii}));
end

% Create a new channels table:
name               = {'V1', 'V2', 'V3', 'V3a', 'V3b', 'hV4', 'LO1','LO2','TO1','IPS'}';
type               = repmat({'n/a'}, [length(name) 1]);
units              = repmat(channels.units(1), [length(name) 1]);
sampling_frequency = repmat(channels.sampling_frequency(1), [length(name) 1]);
subject_name       = subjects;    
number_of_elecs    = nelecs;

channels = table(name, type, units, sampling_frequency, subject_name, number_of_elecs);

end