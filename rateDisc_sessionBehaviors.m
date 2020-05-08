function rateDisc_sessionBehaviors(animal, lastRec)
%% some variables
% nDims = [100, 500, 1000]; %number of whole-frame components to test
% cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %data path to grid server
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server

[dataOverview, ~, ~, ~, ~, ~, ~, ~, trainDates] = rateDiscRecordings; %basic training info
bPath = [cPath animal filesep 'blockData' filesep]; % path for blockdata

sPerformance = []; %nr of local dimensions per session
recDates = [];

%% load explained variance results
load([bPath 'trialInfo.mat'], 'recs');
recIdx = rateDisc_labelRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs);

[Performance,bhv] = rateDisc_sessionPerformance(Animal,cPath,lSessions,highDetection,showPlot,minDelay)


try
    % for session PCs
    load([bPath 'sessionVar.mat'],'sessionVar', 'sessionVarMap', 'sessionAvgVarMap');
    sessionVarMap = arrayShrink(sessionVarMap, allenMask, 'split');
    sessionAvgVarMap = arrayShrink(sessionAvgVarMap, allenMask, 'split');
catch
    errFiles = [];
end



    %for each recording, determine to which part of the training it belongs
    load([bPath 'trialInfo.mat'], 'recs');
    recIdx = rateDisc_labelRecs(trainDates{ismember(dataOverview(:,1), animal)}, recs);
    useIdx = useIdx & recIdx <= lastRec; %only used recordings in assigned range
    useIdx(find(recIdx == lastRec, 1, 'last')) = false; %dont use last session of selected training range
    fprintf('Rejected %d/%d recordings outside of requested training range\n',sum(~(recIdx <= lastRec)),length(useIdx));

    allMse(:,~useIdx) = [];
    sessionVar(:,~useIdx) = [];
    [nDims,b] = sort(nDims);
    allMse = allMse(b,:);
    recs(~useIdx) = [];

    %% get dimensioanlity for local and global dimensions
    for iRecs = 1 : size(allMse,2)
        cData = interp1(nDims,allMse(:,iRecs),nDims(1):nDims(end),'spline');
        gDimCnt(iRecs) = find(cData./max(cData) > dimThresh ,1) + nDims(1) - 1; %number of dimensions need to reach threshold
        lDimCnt(iRecs) = find(sessionVar(:,iRecs)./max(sessionVar(:,iRecs)) > dimThresh ,1); %number of dimensions need to reach threshold
        recDates(iRecs) = datenum(recs(iRecs).name); %date of current recording
    end
    
    
    %% plot overview
    cMap = winter(length(allMse)+1);
    meanMse = mean(allMse,2);
    
    if doPlot
        figure;
        subplot(1,2,1);hold on
        plot(allMse(nDims == 200,:)); axis square; ylabel('mse'); xlabel('recordings'); title(['Prediction error; 200PCs - ' animal])
        xlim([0 sum(useIdx)+1]); ylim([0.95 1]);
        subplot(1,2,2);hold on
        for iFiles = 1 : length(allMse)
            plot(nDims, (allMse(:,iFiles)),'-', 'linewidth',1, 'color',cMap(iFiles+1,:));
            ylabel('explained variance'); xlabel('#PCs'); title(['Dimensionality - ' animal]);
            xlim([0 nDims(end)+nDims(1)]);
        end
        plot(nDims, meanMse,'ko-', 'linewidth',4,'MarkerSize',10); ylim([0.75 1]);
        axis square
        title(['Cumulative explained variance; ' num2str(median(gDimCnt)) ' PCs needed for normalized MSE < ' num2str(1-dimThresh)]);
    end
end
end