function h = SummaryFigure(SessionData,filepath,showData)
%make overview figure and save as pdf

if ~exist('showData', 'var') || isempty(showData)
    showData = false;
end

distFractions = SessionData.DistStim ./ SessionData.TargStim;
distStim = unique(distFractions);
discPerf = [];
UseTrials = SessionData.Assisted;
h = figure(60);

% do general discrimination plot
AxesHandle = subplot(1,2,1); hold on;
for iDist = 1:length(distStim)
    ind1 = distFractions == distStim(iDist) & UseTrials(1:length(distFractions)); %all trials for given distractor
    ind2 = distFractions == distStim(iDist) & SessionData.CorrectSide == 1 & UseTrials(1:length(distFractions)); %left only
    ind3 = distFractions == distStim(iDist) & SessionData.CorrectSide == 2 & UseTrials(1:length(distFractions)); %right only
    
    tCount(iDist,1) = sum(SessionData.Rewarded(ind1)+SessionData.Punished(ind1)); %trialcount
    discPerf(iDist,1) = sum(SessionData.Rewarded(ind1))/tCount(iDist,1); %peformance for given distractor
    
    tCount(iDist,2) = sum(SessionData.Rewarded(ind2)+SessionData.Punished(ind2)); %trialcount
    discPerf(iDist,2) = sum(SessionData.Rewarded(ind2))/tCount(iDist,2); %peformance for given distractor
    
    tCount(iDist,3) = sum(SessionData.Rewarded(ind3)+SessionData.Punished(ind3)); %trialcount
    discPerf(iDist,3) = sum(SessionData.Rewarded(ind3))/tCount(iDist,3); %peformance for given distractor
end

plotDist = distStim(any(~isnan(discPerf),2));
plotPerf = discPerf(any(~isnan(discPerf),2),:);

plot(plotDist,plotPerf(:,2),'-og','Markersize',5,...
    'MarkerEdgeColor','g','MarkerFaceColor','w','linewidth',2,'Parent',AxesHandle)
plot(plotDist,plotPerf(:,3),'-or','Markersize',5,...
    'MarkerEdgeColor','r','MarkerFaceColor','w','linewidth',2,'Parent',AxesHandle)
plot(plotDist,plotPerf(:,1),'-ok','Markersize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','w','linewidth',2,'Parent',AxesHandle)

line([-0.2 1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5], 'parent', AxesHandle)
grid(AxesHandle,'on');
ylim(AxesHandle,[0.45 1.05]); ylabel('Performance');
xlim(AxesHandle,[-0.2 1]); xlabel('Target/Non-target ratio')
legend({'Left','Right','All'},'Location','NorthEast'); axis equal
ylim([0.35 1.05]); title([SessionData.Settings.SubjectName ' - ' date ' - Combined Performance'])

% do modality discrimination plot
discPerf = [];
AxesHandle = subplot(1,2,2); hold on;
ind1 = distFractions == distStim(iDist) & ismember(SessionData.StimType,[3 5:7]) & UseTrials(1:length(distFractions)); %multisensory only
ind2 = distFractions == distStim(iDist) & SessionData.StimType == 1 & UseTrials(1:length(distFractions)); %vision only
ind3 = distFractions == distStim(iDist) & SessionData.StimType == 2 & UseTrials(1:length(distFractions)); %audio only
ind4 = distFractions == distStim(iDist) & SessionData.StimType == 4 & UseTrials(1:length(distFractions)); %audio only

for iDist = 1:length(distStim)
    tCount(iDist,1) = sum(SessionData.Rewarded(ind1)+SessionData.Punished(ind1)); %trialcount
    discPerf(iDist,1) = sum(SessionData.Rewarded(ind1))/tCount(iDist,1); %peformance for given distractor
    
    tCount(iDist,2) = sum(SessionData.Rewarded(ind2)+SessionData.Punished(ind2)); %trialcount
    discPerf(iDist,2) = sum(SessionData.Rewarded(ind2))/tCount(iDist,2); %peformance for given distractor
    
    tCount(iDist,3) = sum(SessionData.Rewarded(ind3)+SessionData.Punished(ind3)); %trialcount
    discPerf(iDist,3) = sum(SessionData.Rewarded(ind3))/tCount(iDist,3); %peformance for given distractor
    
    tCount(iDist,4) = sum(SessionData.Rewarded(ind4)+SessionData.Punished(ind4)); %trialcount
    discPerf(iDist,4) = sum(SessionData.Rewarded(ind4))/tCount(iDist,4); %peformance for given distractor
end

plotDist = distStim(any(~isnan(discPerf),2));
plotPerf = discPerf(any(~isnan(discPerf),2),:);

plot(plotDist,plotPerf(:,1),'-ok','Markersize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','w','linewidth',2,'Parent',AxesHandle)
plot(plotDist,plotPerf(:,2),'-og','Markersize',5,...
    'MarkerEdgeColor','g','MarkerFaceColor','w','linewidth',2,'Parent',AxesHandle)
plot(plotDist,plotPerf(:,3),'-or','Markersize',5,...
    'MarkerEdgeColor','r','MarkerFaceColor','w','linewidth',2,'Parent',AxesHandle)
plot(plotDist,plotPerf(:,4),'-ob','Markersize',5,...
    'MarkerEdgeColor','b','MarkerFaceColor','w','linewidth',2,'Parent',AxesHandle)

line([-0.2 1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5], 'parent', AxesHandle)
grid(AxesHandle,'on');
ylim(AxesHandle,[0.45 1.05]);ylabel('Performance');
xlim(AxesHandle,[-0.2 1]); xlabel('Target/Non-target ratio')
legend({'Mixed','Vision','Audio','Somatosensory'},'Location','NorthEast'); axis equal
ylim([0.35 1.05]); title([SessionData.Settings.SubjectName ' - ' date ' - Modality Performance'])

set(h,'Position',[50 50 1200 800]);
set(h,'PaperOrientation','landscape','PaperPositionMode','auto');
saveas(h,filepath);

if ~showData
    close(h)
end
end