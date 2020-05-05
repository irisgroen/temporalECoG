
conds       = {'CRF', 'ONEPULSE', 'TWOPULSE'};
nChans      = height(channels);

figure;hold on
l = cellfun(@func2str,{results(:).model}, 'UniformOutput', false);

for ii = 1:3
    subplot(1,3,ii); hold on;
    m = [];
    for jj = 1:length(results)
        m(:,jj) = results(jj).R2.concat_cond(ii,:);
        
    end
    h = barh(m);
    set(h, 'BarWidth', 1); 

    numgroups = size(m,1);
    numbars = size(m,2);
    groupwidth = min(0.8,numbars/(numbars+1.5));

    set(gca, 'Xlim', [0 1])
    title(conds{ii});
	set(gca, 'Ylim', [0 nChans+1], 'YTick', 1:nChans, 'YTickLabel', channels.name)% 'YTickLabelRotation', 45);
    ylabel('visual area'); 
    xlabel('R2');
    if ii == 3, legend(l); end
end
set(gcf, 'Position', [400 200 2000 1200]);
