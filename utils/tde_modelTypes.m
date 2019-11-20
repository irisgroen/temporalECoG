function modelTypes = tde_modelTypes()

pth = fullfile(tdeRootPath, 'temporal_models');
d = dir(fullfile(pth, '*.json'));

modelTypes = cell(1,length(d));
for ii = 1:length(d)
    
    [~, fname] = fileparts(d(ii).name);
    modelTypes{ii} = str2func(fname);
end

end