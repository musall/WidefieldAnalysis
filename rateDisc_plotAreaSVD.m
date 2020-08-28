function [gDimCnt, lDimCnt, sessionVar, sessionVarMap, sessionAvgVarMap, modelVarMaps, recs] = rateDisc_plotAreaSVD(cPath, animal, lastRec, doPlot,trainingRange)
%% some variables
% nDims = [100, 500, 1000]; %number of whole-frame components to test
% cPath = '\\grid-hs\churchland_hpc_home\smusall\BpodImager\Animals\'; %data path to grid server
% cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server

bPath = [cPath animal filesep 'blockData' filesep]; % path for blockdata
load([bPath 'mask_' trainingRange '.mat'],'allenMask');

noiseThresh = 0.9; %threshold for mse to reject bad recording
varThresh = 4; %threshold in SDUs for recordings with too high variance (indicative of bad hemo correction)
dimThresh = 0.99; %threshold used to infer dimensionality (default is 1/4 percent reconstruction error below optimum)
gDimCnt = []; %nr of global dimensions per session
lDimCnt = []; %nr of local dimensions per session
recDates = [];

%% load explained variance results
errFiles = dir([bPath 'mse_*' trainingRange '*']);
for iFiles = 1 : length(errFiles)
    temp = textscan(errFiles(iFiles).name, '%s%d%s', 'Delimiter','_');
    nDims(iFiles) = double(temp{2});
    
    load([errFiles(iFiles).folder filesep errFiles(iFiles).name])
    allMse(iFiles,:) = 1-mse;
end

try
    % for session PCs
    load([bPath 'sessionVar_' trainingRange '.mat'],'sessionVar', 'sessionVarMap', 'sessionAvgVarMap');
    sessionVarMap = arrayShrink(sessionVarMap, allenMask, 'split');
    sessionAvgVarMap = arrayShrink(sessionAvgVarMap, allenMask, 'split');
    load([bPath 'modelMaps_' trainingRange '.mat'],'modelVarMaps');
    modelVarMaps = arrayShrink(modelVarMaps.^2, allenMask, 'split');
catch
    errFiles = [];
end

%% check for bad recordings
if ~isempty(errFiles)
    % check for excessive variance or bad reconstruction and don't use recordings
    idx = round(size(sessionVarMap,2) / 2) + (-15 : 15);
    temp = reshape(sessionVarMap(:, idx, :), [], size(sessionVarMap,3)); %get midline activity
    useIdx = allMse(nDims == 200, :) > noiseThresh & zscore(sessionVar(end,:)) < varThresh;

    %for each recording, determine to which part of the training it belongs
    load([bPath 'trialInfo_' trainingRange '.mat'], 'recs');
    for iRecs = 1 : length(recs)
        if isnan(recs(iRecs).name); useIdx(iRecs) = false; end
    end
    fprintf('Rejected %d/%d recordings for mse>%g with 200 PCs or midline variance above %g SDUs\n',sum(~useIdx),length(useIdx),1-noiseThresh, varThresh);
    
    allMse(:,~useIdx) = [];
    sessionVar(:,~useIdx) = [];
    sessionVarMap(:,:,~useIdx) = [];
    sessionAvgVarMap(:,:,~useIdx) = [];
    [nDims,b] = sort(nDims);
    allMse = allMse(b,:);
    recs(~useIdx) = [];
    modelVarMaps(:,:,:,~useIdx) = [];

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