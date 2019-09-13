% old data
load('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess/preproc_BAIRfmri.mat');
Borig = B;
SEorig = SE;

% new data
load('/Volumes/server/Projects/BAIR/Conference/HBM 2019 VisuoTemporal/code/preprocess/preproc_BAIRfmri_10subj.mat');


%% plot single subjects new data
colors = jet(length(subjectList));

figure;hold on
for ii = 1:length(subjectList)
    plot(squeeze(B(ii,1,:))','Color', colors(ii,:),'LineStyle', '-', 'LineWidth', 2);
end
set(gca, 'XTick', 1:17)
set(gca, 'XTickLabel', GLMconditions(trialIndex), 'XTickLabelRotation', 90);
legend(subjectList);
ylabel('beta')
xlabel('stimulus condition');
title('10 subjects')


%% single subjects with SE 
colors = jet(length(subjectList));

figure;
for ii = 1:length(subjectList)
    subplot(2,5,ii); hold on
    %plot(squeeze(B(ii,1,:))','Color', colors(ii,:),'LineStyle', '-', 'LineWidth', 2);
    errorbar(squeeze(B(ii,1,:))',squeeze(SE(ii,1,:))', 'Color', colors(ii,:), 'LineWidth', 2);
    set(gca, 'XTick', 1:17, 'XLim', [1 18], 'YLim', [-1 2.5]);
    set(gca, 'XTickLabel', GLMconditions(trialIndex), 'XTickLabelRotation', 90);

    title(subjectList(ii))
    if ismember(ii, [2 3 4])
        errorbar(squeeze(Borig(ii-1,1,:))',squeeze(SEorig(ii-1,1,:))', 'Color', colors(ii,:), 'LineStyle', ':', 'LineWidth', 2);
        legend({'new', 'original'});
    end
end

%% compare 

Bnew = B(2:4,:,:);
SEnew = SE(2:4,:,:);
subjectListNew = subjectList(2:4);
colors = [1 0 0; 0 1 0; 0 0 1];

figure; hold on
c = 1; names = [];
for ii = 1:3
    plot(squeeze(Borig(ii,1,:))', 'Color', colors(ii,:), 'LineStyle', '--', 'LineWidth', 2);
    names{c} = [subjectListNew{ii} ' original']; c = c+1;
    plot(squeeze(Bnew(ii,1,:))', 'Color', colors(ii,:), 'LineStyle', '-', 'LineWidth', 2);
	names{c} = [subjectListNew{ii} ' new']; c = c+1;
end
set(gca, 'XTick', 1:17)
set(gca, 'XTickLabel', GLMconditions(trialIndex), 'XTickLabelRotation', 90);
legend(names);
ylabel('beta')
xlabel('stimulus condition');
