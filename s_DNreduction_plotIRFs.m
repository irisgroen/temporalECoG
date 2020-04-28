figure;hold on;
subplot(1,3,1);hold on;
plot(t,irf, 'k', 'LineWidth', 2)
plot(t,irf_pos, 'b', 'LineWidth', 2)
plot(t,irf_neg, 'c', 'LineWidth', 2)
plot(t,-prm.weight.*irf_neg, 'c:', 'LineWidth', 2)
title('DN')
xlim([0 1])

subplot(1,3,2);hold on;
plot(t,irf, 'k', 'LineWidth', 2)
plot(t,irf_pos, 'b', 'LineWidth', 2)
plot(t,irf_neg, 'c', 'LineWidth', 2)
plot(t,-prm.weight.*irf_neg, 'c:', 'LineWidth', 2)
title('LINEAR_RECTF_EXP_NORM_DELAY')
xlim([0 1])

subplot(1,3,3);hold on;
plot(t,irf, 'k', 'LineWidth', 2)
plot(t,irf_pos, 'b', 'LineWidth', 2)
plot(t,irf_neg, 'c', 'LineWidth', 2)
plot(t,-prm.weight.*irf_neg, 'c:', 'LineWidth', 2)
title('LINEAR_RECTH_EXP_NORM_DELAY')
xlim([0 1])

%% irf-pos
nROI = size(results(1).params,2);
t = 0:0.001:3;

figure;hold on;
for ii = 1:nROI
    l = [];
    subplot(3, nROI/3, ii);hold on;
    %DN
    tau = results(1).params(1,ii);
    n = 2;
    y = gammaPDF(t,tau,n);
    plot(t,y,'b', 'LineWidth', 2);
    l{1} = [tau n];
    % RECTF
    tau = results(2).params(1,ii);
    n = max(round(results(2).params(2,ii)),1);
    y = gammaPDF(t,tau,n);
    plot(t,y,'r', 'LineWidth', 2);
    l{2} = [tau n];
    % RECTH
    tau = results(3).params(1,ii);
    n = max(round(results(3).params(2,ii)),1);
    y = gammaPDF(t,tau,n);
    plot(t,y,'g', 'LineWidth', 2);
    l{3} = [tau n];
    legend({sprintf('tau1 = %0.2f n1 = %d', l{1}(1),l{1}(2)), ...
            sprintf('tau1 = %0.2f n1 = %d', l{2}(1),l{2}(2)), ...
            sprintf('tau1 = %0.2f n1 = %d', l{3}(1),l{3}(2))});
	title(channels.name{ii});
    xlim([0 0.5]);
    ylim([0 0.02]);
end
set(gcf, 'Position', [ 25 25 1250 750]);

%% irf-neg
nROI = size(results(1).params,2);
t = 0:0.001:3;

figure;hold on;
for ii = 1:nROI
    l = [];
    subplot(3, nROI/3, ii);hold on;
    %DN
    tau = results(1).params(1,ii)*1.5;
    n = 2;
    w = results(1).params(2,ii);
    y = -w .* gammaPDF(t,tau,n);
    plot(t,y,'b', 'LineWidth', 2);
    l{1} = [tau n w];
    % RECTF
    tau = results(2).params(3,ii);
    n = max(round(results(2).params(4,ii)),1);
	w = results(2).params(5,ii);
    y = -w .* gammaPDF(t,tau,n);
    plot(t,y,'r', 'LineWidth', 2);
    l{2} = [tau n w];
    % RECTH
    tau = results(3).params(3,ii);
    n = max(round(results(3).params(4,ii)),1);
	w = results(3).params(5,ii);
    y = -w .* gammaPDF(t,tau,n);
    plot(t,y,'g', 'LineWidth', 2);
    l{3} = [tau n w];
    legend({sprintf('tau2 = %0.2f n2 = %d w = %0.2f', l{1}(1),l{1}(2), l{1}(3)), ...
            sprintf('tau2 = %0.2f n2 = %d w = %0.2f', l{2}(1),l{2}(2), l{2}(3)), ...
            sprintf('tau2 = %0.2f n2 = %d w = %0.2f', l{3}(1),l{3}(2), l{3}(3))});
	title(channels.name{ii});
    xlim([0 0.5]);
	ylim([-0.005 0.005]);

end
set(gcf, 'Position', [ 25 25 1250 750]);

%% irf
nROI = size(results(1).params,2);
t = 0:0.001:3;

figure;hold on;
for ii = 1:nROI
    l = [];
    subplot(3, nROI/3, ii);hold on;
    %DN
    tau1 = results(1).params(1,ii);
    n1 = 2;
    w = results(1).params(2,ii);
    irf_pos = gammaPDF(t, tau1, n1);
    irf_neg = gammaPDF(t, tau1*1.5, n1);
    y = irf_pos - w.* irf_neg;
    plot(t,y,'b', 'LineWidth', 2);
    % RECTF
    tau1 = results(2).params(1,ii);
    tau2 = results(2).params(3,ii);
    n1 = max(round(results(2).params(2,ii)),1);
	n2 = max(round(results(2).params(4,ii)),1);
	w = results(2).params(5,ii);
    irf_pos = gammaPDF(t, tau1, n1);
    irf_neg = gammaPDF(t, tau2, n2);
    y = irf_pos - w.* irf_neg;
    plot(t,y,'r', 'LineWidth', 2);
    % RECTH
    tau1 = results(3).params(1,ii);
    tau2 = results(3).params(3,ii);
    n1 = max(round(results(3).params(2,ii)),1);
	n2 = max(round(results(3).params(4,ii)),1);
	w = results(3).params(5,ii);
    irf_pos = gammaPDF(t, tau1, n1);
    irf_neg = gammaPDF(t, tau2, n2);
    y = irf_pos - w.* irf_neg;
    plot(t,y,'g', 'LineWidth', 2);
    if ii==1, legend('DN', 'LINEAR_RECTF_EXP_NORM_DELAY', 'LINEAR_RECTH_EXP_NORM_DELAY');end
	title(channels.name{ii});
    xlim([0 1]);
    ylim([-0.005 0.02]);

end
set(gcf, 'Position', [ 25 25 1250 750]);