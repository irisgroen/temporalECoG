
% Load data and fits
modelfun = @LINEAR_RECTF_EXP_NORM_DELAY;%@DN;%
xvalmode = 0;
datatype = 'individualelecs'; % 'electrodeaverages';

datastr = 'bads';
[D] = tde_loadDataForFigure(modelfun, xvalmode, datatype, [], datastr);
[results_bads] = tde_evaluateModelFit(D,0);

% datastr = 'surrogateopts';
% [D] = tde_loadDataForFigure(modelfun, xvalmode, datatype, [], datastr);
% [results_surr] = tde_evaluateModelFit(D,0);

% datastr = 'lsqnonlin';
% [D] = tde_loadDataForFigure(modelfun, xvalmode, datatype, [], datastr);
% [results_lsq] = tde_evaluateModelFit(D,0);
% 
% datastr = 'fmincon';
% [D] = tde_loadDataForFigure(modelfun, xvalmode, datatype, [], datastr);
% [results_fmincon] = tde_evaluateModelFit(D,0);

p = loadjson(fullfile(tdeRootPath, 'temporal_models', sprintf('%s.json', func2str(results_bads.model))));
T = strsplit(p.params, ',');

figure, hold on,
subplot(2,5,1); hold on
%plot(results_bads.R2.concat_all,'-o', 'LineWidth', 2);
%plot(results_surr.R2.concat_all, '-.s','LineWidth', 2);
%plot(results_lsq.R2.concat_all, '--d', 'LineWidth', 2);
%plot(results_fmincon.R2.concat_all, ':x', 'LineWidth', 2);
plot(results_bads.R2.concat_all,'LineWidth', 2);
%plot(results_lsq.R2.concat_all, 'LineWidth', 2);
%plot(results_fmincon.R2.concat_all,   'LineWidth', 2);

ylim([0 1]);

%legend('bads','surrogateopts', 'lsqnonlin', 'fmincon'); xlabel('electrode'); ylabel('fitted value');
legend('bads','lsqnonlin', 'fmincon'); xlabel('electrode'); ylabel('fitted value');

title('R2');

for ii = 1:length(T)
    subplot(2,5,ii+1); hold on;
    %plot(results_bads.params(ii,:),'-o', 'LineWidth', 2);
	%plot(results_surr.params(ii,:), '--d', 'LineWidth', 2);
    %plot(results_lsq.params(ii,:), '-.s', 'LineWidth', 2);
    %plot(results_fmincon.params(ii,:),  ':x','LineWidth', 2);
    plot(results_bads.params(ii,:), 'LineWidth', 2);
   % plot(results_lsq.params(ii,:), 'LineWidth', 2);
   % plot(results_fmincon.params(ii,:), 'LineWidth', 2);
    title(T{ii});
end

set(findall(gcf,'-property','FontSize'),'FontSize',18)
set(gcf, 'position',  get(0, 'screensize'));
