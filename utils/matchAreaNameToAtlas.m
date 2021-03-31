function area_idx = matchAreaNameToAtlas(areaName, atlasLabels)
% find atlasLabels for that contain a given areaName (e.g. V1v and V1d
% for V1 for wang atlases), preventing matches for V3a/b with V3 and
% collapsing all IPS maps


if strcmpi(areaName, 'V3')
    area_idx = contains(atlasLabels, 'V3') & ~contains(atlasLabels, {'V3a', 'V3b'});
elseif strcmpi(areaName, 'IPS')
    area_idx = contains(atlasLabels, {'IPS0', 'IPS1', 'IPS2', 'IPS3', 'IPS4', 'IPS5'});
elseif strcmpi(areaName, 'V3ab')
    area_idx = contains(atlasLabels, {'V3a', 'V3b'});
elseif strcmpi(areaName, 'LOTO')
    area_idx = contains(atlasLabels, {'LO1', 'LO2', 'TO1', 'TO2'});
elseif strcmpi(areaName, 'TO')
    area_idx = contains(atlasLabels, {'TO1', 'TO2'});
elseif strcmpi(areaName, 'V123')
    area_idx = contains(atlasLabels, {'V1', 'V2', 'V3'}) & ~contains(atlasLabels, {'V3a', 'V3b'});
elseif strcmpi(areaName, 'higher')
    area_idx = contains(atlasLabels, {'V3a', 'V3b','LO1', 'LO2', 'TO1', 'TO2', 'IPS0', 'IPS1', 'IPS2', 'IPS3', 'IPS4', 'IPS5'});
else
    area_idx = contains(atlasLabels, areaName); 
end

end