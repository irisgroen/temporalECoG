
conditionsOfInterest = [5 11];
normalize = 0;

% Plot multiple areas together in one plot
nChans = height(channels);
colors = jet(nChans);

subjectNames = unique(channels.subject_name);
%subjectNames = {'chaam', 'som692', 'som708'};
nSub = length(subjectNames);

if normalize 
    figName = 'FullContrast500ms normalized';
else
    figName = 'FullContrast500ms';
end
figure('Name', figName);hold on    

for kk = 1:nSub
    
    % Determine how many subplots to make
    nRow  = ceil(sqrt(nSub));
    nCol  = ceil(sqrt(nSub));
    if nSub<= (nRow*nCol)-nCol, nRow = nRow-1;end

    subIdx = find(contains(channels.subject_name, subjectNames{kk}));
    subdata = data2fit(:,:,subIdx);
    subchans = channels(subIdx,:);
    %[INX,~] = groupElecsByVisualArea(subchans);
    [subdata, subdatase, subchans] = average_elecs(subdata, subchans, conditionsOfInterest, 0);
    nChans = height(subchans);
    colors = jet(nChans);
    
    subplot(nRow,nCol,kk); hold on
    %subplot(1,3,kk); hold on
    for ii = 1:height(subchans)   
        %m = squeeze(mean(subdata(:,conditionsOfInterest,ii),2));
        m = subdata(:,ii);
        %sm = squeeze(mean(subdata(:,conditionsOfInterest,INX{ii}),2));
        if normalize, m = m./max(m); end
        if ~isempty(m)
            %plot(t,m,'Color', colors(ii,:), 'LineWidth', 2);        
            if normalize
                ecog_plotSingleTimeCourse(t, m, [], colors(ii,:), [], 'Response normalized', [-0.1 1.1], 'Time');
            else
                ecog_plotSingleTimeCourse(t, m, [], colors(ii,:), [], 'Response', [], 'Time');
            end
            %ecog_plotSingleTimeCourse(t, m, [], colors(ii,:));

        end
        
        %se = squeeze(mean(subdatase(:,conditionsOfInterest,ii),2));
        %h = ciplot(m(kk,:,ii)-se(kk,:,ii), m(kk,:,ii)+se(kk,:,ii), [], colors{kk}, 0.25);
        %h = ciplot(se(kk,:,ii,1), se(kk,:,ii,2), [], colors{kk}, 0.25);
        %h.Annotation.LegendInformation.IconDisplayStyle = 'off';      
    end
    legend(subchans.name); 
    %xlabel('Time'), title('Data'); ylabel('Response');
    %if normalize, set(gca, 'Ylim', [-0.1 1.1]);ylabel('Response (normalized)');end
    title(sprintf('FullContrast500ms %s', subjectNames{kk}));
    
%     %subplot(1,3,2); hold on
%     for ii = 1:nChans    
%         m = squeeze(mean(results(1).pred(:,conditionsOfInterest,ii),2));
%         if normalize, m = m/max(m); end
%         plot(t,m,'Color', colors(ii,:), 'LineWidth', 2);        
%     end
%     legend(channels.name); xlabel('Time'), title('DN Prediction');
% 
%     if normalize, set(gca, 'Ylim', [-0.1 1.1]);end
% 
%     subplot(1,3,3); hold on
%     for ii = 1:nChans    
%         m = results(1).derived.pred(:,ii);
%         if normalize, m = m/max(m); end
%         plot(m,'Color', colors(ii,:), 'LineWidth', 2);        
%     end
%     xlabel('Time'), title('DN derivedPrediction');
%     %set(gca, 'FontSize', 14);
% 
%     set(gca, 'Xlim', [0 1000]);
%     if normalize, set(gca, 'Ylim', [-0.1 1.1]);end
% 
%     set(gcf, 'Position', [600 700 1600 600]);
 set(gca, 'YLim', [0 40]);
end
set(findall(gcf,'-property','FontSize'),'FontSize',14)
set(gcf, 'Position', [400 200 2000 1200]);

% Data averaging
function [mdata, sedata, channels] = average_elecs(data, channels, conditionsOfInterest, normalize)
    
    [INX, channels] = groupElecsByVisualArea(channels);
    
    for ii = 1:length(INX)
        chanInx = find(INX{ii});
        data_to_average = squeeze(mean(data(:,conditionsOfInterest,chanInx),2)); % average conditions first
        if normalize 
            for jj = 1:length(chanInx)
                data_to_average(:,jj) = data_to_average(:,jj) ./ max(data_to_average(:,jj));
            end 
        end
        m = mean(data_to_average,2);
        se = std(data_to_average,0,2)/sqrt(length(chanInx));
        mdata(:,ii) = m;
        sedata(:,ii) = se;
    end
end