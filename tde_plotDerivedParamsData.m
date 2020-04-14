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
stim_on = t>0 & t<0.5;

figure('Position', [360    44   879   654]); hold on;

for ii = 1:length(conditionsOfInterest)
    subplot(2,2,ii); hold on
	inx = find(contains(stim_info.name, conditionsOfInterest{ii}));
    m = squeeze(sum(data(stim_on, inx, :),1));
    m = m./max(m);
	cmap = num2cell(parula(size(m,2)),2);
    h = plot(m, '.-', 'MarkerSize', 20, 'LineWidth', 2);
    xlim([0 length(inx)+1]);
    ylim([0 1.1]);
    set(gca, 'XTick', 1:length(inx));
    if ii == 1, set(gca, 'XTickLabel', stim_info.contrast(inx)); xlabel('contrast level'); title('contrast response functions'); end
    if ii == 2, set(gca, 'XTickLabel', stim_info.duration(inx)); xlabel('duration (s)'); title('temporal summation onepulse'); end
    if ii == 3, set(gca, 'XTickLabel', stim_info.ISI(inx)); xlabel('ISI (s)'); title('temporal summation twopulse'); end
	set(h, {'color'}, cmap);
    ylabel('summed broadband (0-0.5s)');
end

[m] = tde_computeISIrecovery(data,channels,t,stim_info);
stimISI = stim_info.ISI(inx);

subplot(2,2,ii+1); hold on
cmap = num2cell(parula(size(m,2)),2);
h = plot(stimISI, m, '.-', 'MarkerSize', 20, 'LineWidth', 2);
xlabel('ISI (s)');
ylabel('percentage of first pulse');
set(h, {'color'}, cmap);
xlim([stimISI(1)-0.1 stimISI(end)+0.1]);
ylim([0 1]);
legend(channels.name, 'Location', 'SouthEast');
title('recovery with ISI - all areas');     

end
