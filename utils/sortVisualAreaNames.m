function [names_sorted,sort_index] = sortVisualAreaNames(names)


area_labels = {'V1', 'V1v', 'V1d', ...
               'V2', 'V2v', 'V2d', ...
               'V3', 'V3v', 'V3d', ...
               'V3a', 'V3b', ...
               'hV4', ...
               'LO1', 'LO2', ...
               'TO1', 'TO2', ...
               'VO1', 'VO2', ...
               'IPS0', 'IPS1', 'IPS2', 'IPS3', 'IPS4', 'IPS5' ...
               'PHC1', 'PHC2', ... 
               'SPL1','FEF', ...
               'none'}; 

sortorder = nan(size(names));

for ii = 1:length(names)
    inx = strcmp(names{ii},area_labels);
    %inx = contains(area_labels, names{ii});
    if max(inx) > 0
        matchedArea = area_labels{inx};
        sortorder(ii) = find(strcmp(matchedArea,area_labels));
    end
end

[~,sort_index] = sort(sortorder);
names_sorted = names(sort_index);

end