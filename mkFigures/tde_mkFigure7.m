% tde_mkFigure 7


modelfuns = @LINEAR_RECTF_EXP_NORM_DELAY;
xvalmode = 0;
datatype = 'individualelecs';

for ii = 1:length(modelfuns)
    [d(ii)] = tde_loadDataForFigure(modelfuns{ii}, xvalmode, datatype);
end

[results] = tde_evaluateModelFit(d);