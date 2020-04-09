function tde_plotDerivedParamsData(data,channels,t,stim_info)

% For each channel in the data, compute 3 things:
% - contrast response function (based on CRF trials)
% - subadditive temporal responses (based on ONEPULSE trials)
% - subadditive temporal responses (based on TWOPULSE trials)
% - time to recovery (based on TWOPULSE trials)



[nSamp, nChan, nDatasets] = size(data);

% Determine if data was averaged across elecs prior to fit
if isfield(summary(channels), 'number_of_elecs')
    dataWasAveraged = true;
else
    dataWasAveraged = false;
end

conditionsOfInterest = {'CRF', 'ONEPULSE', 'TWOPULSE'};
stim_on = t>0 && t<0.5;

for ii = 1:length(conditionsOfInterest)
    subplot(2,2,ii); hold on
	inx = contains(stim_info.name, conditionsOfInterest{jj});
    m = squeeze(sum(data(stim_on, inx, :),1));
	cmap = num2cell(parula(size(m,2)),2);
    h = plot(1:length(inx), m, '.-', 'MarkerSize', 20, 'LineWidth', 2);
	set(h, {'color'}, cmap);
end

subplot(2,2,ii+1); hold on
    
    
        

end
