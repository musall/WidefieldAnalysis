%% basic variables
cPath = 'U:\space_managed_data\BpodImager\Animals\'; %data path on the server
trialSegments = [0 55 80 130 160]; %segments to create maps for whole trial regressors
motorIdx = 16; %index for zero-lag motor regressor
dimCnt = 200; % use only dimCnt dimensions from Vc. This should match what was used in the model
baseRange = 1:15; %baseline frames
areaIdx = {'VISp' 'RSPd' 'SSp-ll' 'MOs'}; % V1 - RS - HL - M2
normData = false; %flag to zscore data

%% session data
[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, segIdxRealign] = delayDecRecordings;
taskLabels = [sensorLabels cogLabels]; %combine sensory and cognitive labels, so sensory is really 'non-motor' or 'task' here
segIdxRealign{2} = 46:65;
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts
cSegs = cat(2,segIdxRealign{ismember(segLabels,{'Stim1' 'Stim2'})}); %time index for stimulus segments
cSegs = [cSegs(1) segIdxRealign{ismember(segLabels,{'Wait'})}(round(end/2))]; %add delay period
cSegs = cSegs(1) : cSegs(end);

%% get reconstructed Vs, modality index, datapath and allOpts
[recV, recLabels, dataPath, allOpts, modIdx, sideIdx, alignIdx, baseLength, frames, stimTimes] = BpodImager_motorReconstruct('All', 'All', false); %get reconstructed V, used full model
recV = cat(5,recV{:});
stimTimes = cat(2,stimTimes{:});
stimTimes(stimTimes > 1.5) = [];
figure; histogram(stimTimes,50); axis square; %show stim time histogram
nRecLabels = [{'full'} recLabels{:} {'motor' 'motorVideo' 'sensory' 'cognitive'}];
altRecLabels = nRecLabels; %same as nReclabels but with expReward and novReward
altRecLabels(ismember(nRecLabels,'visReward')) = {'expReward'};
altRecLabels(ismember(nRecLabels,'audReward')) = {'novReward'};
                
 %get allen maps
load('allenDorsalMap.mat')
addEdgeOutlinesToDorsalMap(dorsalMaps); %make area figure for orientation
mask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
[x1, y1] = size(mask);
load([dataPath{1} 'snapshot_1.mat'])
snap = alignAllenTransIm(single(snap),allOpts(1).opts.transParams);
[x2, y2] = size(snap);
mask = mask(1:min([x1 x2]), 1:min([y1 y2])); %cut mask to size
rightHs = ismember(dorsalMaps.sidesSplit,'R'); %index for labels on the left HS

% get allen area masks
for x = 1:length(areaIdx)
    ind = ismember(dorsalMaps.labelsSplit,areaIdx{x}) & rightHs;
    areaCoord{x} = poly2mask(dorsalMaps.edgeOutlineSplit{ind}(:,2), dorsalMaps.edgeOutlineSplit{ind}(:,1),size(mask,1),size(mask,2));
end
rightHs = find(rightHs); % merge outlines
meanDev = NaN(length(taskLabels)-1, length(dataPath), 'single');

%% load some raw data and compute psth deviations
tic;
for iAnimals = 1:length(dataPath)
    load([dataPath{iAnimals} 'Vc.mat'],'U')
    U = U(:,:,1:dimCnt);
    U = alignAllenTransIm(U,allOpts(iAnimals).opts.transParams);
    U = arrayShrink(U(1:size(mask,1),1:size(mask,2),:),mask);
    
    load([dataPath{iAnimals} 'dimBeta.mat'],'dimBeta');
    load([dataPath{iAnimals} 'regData.mat'],'recIdx', 'idx');
    cIdx = recIdx(~idx);
    
    lVis = dimBeta(ismember(cIdx, find(ismember(recLabels, 'lVisStim'))), :);
    rVis = dimBeta(ismember(cIdx, find(ismember(recLabels, 'rVisStim'))), :);
    cInd = 1 : min([size(lVis,1) size(rVis,1)]);
    vis = lVis(cInd, :) + rVis(cInd, :); clear lVis rVis
    
    
    lAud = dimBeta(ismember(cIdx, find(ismember(recLabels, 'lAudStim'))), :);
    rAud = dimBeta(ismember(cIdx, find(ismember(recLabels, 'rAudStim'))), :);
    cInd = 1 : min([size(lAud,1) size(rAud,1)]);
    aud = lAud(cInd, :) + rAud(cInd, :); clear lAud rAud
    
    choice = dimBeta(ismember(cIdx, find(ismember(recLabels, 'Choice'))), :);
    time = dimBeta(ismember(cIdx, find(ismember(recLabels, 'time'))), :);
    if visExp(iAnimals)
        reward = dimBeta(ismember(cIdx, find(ismember(recLabels, 'visReward'))), :);
    else
        reward = dimBeta(ismember(cIdx, find(ismember(recLabels, 'audReward'))), :);
    end
    
    if iAnimals == 1
        visBeta = NaN(sum(~mask(:)), 120, length(dataPath), 'single');
        audBeta = NaN(sum(~mask(:)), 120, length(dataPath), 'single');
        choiceBeta = NaN(sum(~mask(:)), 200, length(dataPath), 'single');
        rewardBeta = NaN(sum(~mask(:)), 200, length(dataPath), 'single');
        timeBeta = NaN(sum(~mask(:)), 120, length(dataPath), 'single');
    end
    
    visBeta(:, 1 : size(vis,1), iAnimals) = U * vis';
    audBeta(:, 1 : size(aud,1), iAnimals) = U * aud';
    choiceBeta(:, 1 : size(choice,1), iAnimals) = U * choice';
    rewardBeta(:, 1 : size(reward,1), iAnimals) = U * reward';
    timeBeta(:, 1 : size(time,1), iAnimals) = U * time';

    % compute mean overall deviation
    for iRegs = 1 : length(taskLabels) - 1
        temp = dimBeta(ismember(cIdx, find(ismember(recLabels, taskLabels(iRegs)))), :);
        if any(ismember({'lVisStim' 'rVisStim' 'lAudStim' 'rAudStim'},taskLabels{iRegs}))
            temp = arrayShrink(U * temp(1 : length(cSegs),:)', mask, 'split');
        else
            temp = arrayShrink(U * temp(cSegs,:)', mask, 'split');
        end
        
        %check for expertise for reward regressors
        if (visExp(iAnimals) && strcmpi(taskLabels(iRegs), 'visReward')) || (audExp(iAnimals) && strcmpi(taskLabels(iRegs), 'audReward'))
            meanDev(strcmpi(taskLabels, 'visReward'), iAnimals) = nanmean(abs(temp(:)));
        elseif (visExp(iAnimals) && strcmpi(taskLabels(iRegs), 'audReward')) || (audExp(iAnimals) && strcmpi(taskLabels(iRegs), 'visReward'))
            meanDev(strcmpi(taskLabels, 'audReward'), iAnimals) = nanmean(abs(temp(:)));
        else
            meanDev(iRegs, iAnimals) = nanmean(abs(temp(:)));
        end
    end
    
    if rem(iAnimals, round(length(dataPath)/5)) == 0
        fprintf(1, 'Done. Current recording is %d out of %d\n', iAnimals, length(dataPath));
        toc;
    end
end


%% single regressor contributions
idx = ones(1,length(taskLabels)-1);

figure
ax = regressorPlot(meanDev(idx>0,:)',taskLabels(idx>0),5,subplot(2,2,1:2),[0 1 0],idx(idx>0),0.6);


%% compute psth deviations
stimAlign{1} = 1 : 18;
stimAlign{2} = 33 : 51;
stimAlign{3} = 52 : 65;

figure
subplot(2,3,1)
cMovie = nanmean(arrayShrink(abs(visBeta),mask,'split'),4);
cRange = [0 prctile(cMovie(:),85)];
mapImg = imshow(nanmean(cMovie,3), cRange);
colormap(mapImg.Parent,jet(256)); axis image
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
title('Vision')

subplot(2,3,2)
cMovie = nanmean(arrayShrink(abs(audBeta),mask,'split'),4);
% cRange = [0 prctile(cMovie(:),90)];
mapImg = imshow(nanmean(cMovie,3), cRange);
colormap(mapImg.Parent,jet(256)); axis image
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
title('Audio')

subplot(2,3,3)
cMovie = nanmean(arrayShrink(abs(timeBeta),mask,'split'),4);
% cRange = [0 prctile(cMovie(:),95)];
mapImg = imshow(nanmean(cMovie,3), cRange);
colormap(mapImg.Parent,jet(256)); axis image
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
title('Time')

subplot(2,3,4)
cMovie = nanmean(arrayShrink(abs(rewardBeta(:, cSegs, :)),mask,'split'),4);
cRange = [0 prctile(cMovie(:),90)];
mapImg = imshow(nanmean(cMovie,3), cRange);
colormap(mapImg.Parent,jet(256)); axis image
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
title('Reward')

subplot(2,3,5)
cMovie = nanmean(arrayShrink(abs(choiceBeta(:, cSegs, :)),mask,'split'),4);
% cRange = [0 prctile(cMovie(:),90)];
mapImg = imshow(nanmean(cMovie,3), cRange);
colormap(mapImg.Parent,jet(256)); axis image
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
title('Choice')




%    %% 
% figure;
% for iSegs = 1:3
%     subplot(1,3,iSegs)
%     cMap = abs(nanmean(cMovie(:,:,stimAlign{iSegs}),3)/2);
%     cRange = [0 prctile(cMap(:),90)];
% 
%     mapImg = imshow(cMap, cRange);
%     colormap(mapImg.Parent,jet(256)); axis image
%     set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
%     title(['All mice. cVision: ' segLabels{iSegs+1}])
%     
%     hold(mapImg.Parent, 'on');
%     for x = 1: length(rightHs)
%         plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
%     end
% end




