weksel%% basic variables
cPath = 'U:\space_managed_data\BpodImager\Animals\'; %data path on the server
trialSegments = [0 55 80 130 160]; %segments to create maps for whole trial regressors
motorIdx = 16; %index for zero-lag motor regressor
dimCnt = 200; % use only dimCnt dimensions from Vc. This should match what was used in the model
baseRange = 1:15; %baseline frames
areaIdx = {'VISp' 'RSPd' 'SSp-ll' 'MOs'}; % V1 - RS - HL - M2
normData = false; %flag to zscore data
iMod = 1; % 3 = all trials
iSeg = 4; % 4 = all frames

%% session data
[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, segIdxRealign] = delayDecRecordings;
sensorLabels = [sensorLabels cogLabels]; %combine sensory and cognitive labels, so sensory is really 'non-motor' or 'task' here
segIdxRealign{2} = 46:65;
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts

%% get reconstructed Vs, modality index, datapath and allOpts
[recV, recLabels, dataPath, allOpts, modIdx, sideIdx, alignIdx, baseLength, frames, stimTimes] = BpodImager_motorReconstruct('All', 'All', false); %get reconstructed V, used full model
recV = cat(5,recV{:});

%load extra group labels
load([dataPath{1} 'predVariance' filesep 'extraGroups.mat'],'extraGroups');
extraGroups{1,end + 1} = 'task';
extraGroups{2,end} = recLabels(~ismember(recLabels,motorLabels));
extraGroups{1,end + 1} = 'full';
extraGroups{2,end} = recLabels;

%load other motor labels
load([dataPath{1} 'predVariance' filesep 'oMotorLabels.mat'],'oMotorLabels'); 
fullRecLabels = [recLabels extraGroups(1,:)];

% stimTimes = cat(2,stimTimes{:});
% stimTimes(stimTimes > 1.5) = [];
% figure; histogram(stimTimes,50); axis square; %show stim time histogram

%rearrange recV for expert and novice reward (expReward = visReward, audReward = novReward)
temp = recV;
recV(:,:,:,audExp,ismember(recLabels,'audReward')) = temp(:,:,:,audExp,ismember(recLabels,'visReward'));
recV(:,:,:,audExp,ismember(recLabels,'visReward')) = temp(:,:,:,audExp,ismember(recLabels,'audReward'));

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
xRange = [(-baseLength : 1 : 0)./30, (1 : frames-baseLength-1) ./30]; %time vector for x-axis. Stimulus onset is at 0s.

%% get reconstructed data
fullRec{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
motorRec{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
motorRec{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
sensoryRec{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
sensoryRec{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
allData{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

% frame index for stimulus and waiting segments
cSegs = cat(2,segIdxRealign{ismember(segLabels,{'Stim1' 'Stim2'})}); %time index for stimulus segments
cSegs = [cSegs(1) segIdxRealign{ismember(segLabels,{'Wait'})}(round(end/2))]; %add delay period
cSegs = cSegs(1):cSegs(end);
respVar = NaN(sum(~mask(:)), 3, length(fullRecLabels), length(dataPath), 4, 'single');

% convolve real and reconstructed Vs to create PSTHs
tic
for iAnimals = 1:length(dataPath)
    
    load([dataPath{iAnimals} 'Vc.mat'],'U')
    U = U(:,:,1:dimCnt);
    U = alignAllenTransIm(U,allOpts(iAnimals).opts.transParams);
    U = arrayShrink(U(1:size(mask,1),1:size(mask,2),:),mask);
    load([dataPath{iAnimals} 'interpVc.mat'],'Vc')
    Vc = bsxfun(@minus, Vc, mean(Vc,2)); %make sure Vc is zero-mean

    Vc = Vc(:,alignIdx{iAnimals});
    Vc = mean(reshape(Vc(:,modIdx{iAnimals,iMod}),dimCnt,frames,[]),3);
    
    allData{1}(:, :, iAnimals) = U * Vc; % original data
    fullRec{1}(:, :, iAnimals) = U * sum(recV(:,:,iMod,iAnimals,:),5); % full reconstruction

    motorRec{1}(:, :, iAnimals) = U * sum(recV(:,:,1,iAnimals,ismember(recLabels,motorLabels)),5); % motor reconstruction
    motorRec{2}(:, :, iAnimals) = U * sum(recV(:,:,2,iAnimals,ismember(recLabels,motorLabels)),5); % motor reconstruction
    motorRec{3}(:, :, iAnimals) = U * sum(recV(:,:,iMod,iAnimals,ismember(recLabels,motorLabels)),5); % motor reconstruction

    sensoryRec{1}(:, :, iAnimals) = U * sum(recV(:,:,1,iAnimals,ismember(recLabels,sensorLabels)),5); % non-motor reconstruction
    sensoryRec{2}(:, :, iAnimals) = U * sum(recV(:,:,2,iAnimals,ismember(recLabels,sensorLabels)),5); % non-motor reconstruction
    sensoryRec{3}(:, :, iAnimals) = U * sum(recV(:,:,iMod,iAnimals,ismember(recLabels,sensorLabels)),5); % motor reconstruction

    % Check for ramp angle during stimulus and waiting period in each regressor
    for iMods = 1:3
        for iRegs = 1 : length(fullRecLabels)
                
            if any(strcmp(recLabels, fullRecLabels{iRegs}))
                cIdx = find(strcmp(recLabels, fullRecLabels{iRegs}));
            else
                cIdx = ismember(recLabels, extraGroups{2, ismember(extraGroups(1,:), fullRecLabels{iRegs})});
            end
            
            %determine trained/untrained modality for iMods 1/2. 3 is full PSTH.
            if iMods == 1 && visExp(iAnimals) || iMods == 2 && ~visExp(iAnimals)
                cMovie = U * sum(recV(:,:,1,iAnimals,cIdx),5);
            elseif iMods == 1 && ~visExp(iAnimals) || iMods == 2 && visExp(iAnimals)
                cMovie = U * sum(recV(:,:,2,iAnimals,cIdx),5);
            elseif iMods == 3
                cMovie = U * sum(recV(:,:,3,iAnimals,cIdx),5);
            end
            
            for iSegs = 1:3
                respVar(:,iMods,iRegs,iAnimals,iSegs) = mean(abs(cMovie(:,segIdxRealign{iSegs+2})),2); %signal deviation from zero during stimulus/delay period
            end
            respVar(:,iMods,iRegs,iAnimals,iSegs+1) = mean(abs(cMovie),2); %signal deviation from zero during stimulus/delay period
            respVar(:,iMods,iRegs,iAnimals,iSegs+2) = mean(abs(cMovie(:,cSegs)),2); %signal deviation from zero during stimulus/delay period
            
            clear cMovie cBetas
        end
    end
    
    if rem(iAnimals, round(length(dataPath)/5)) == 0
        fprintf(1, 'Done. Current recording is %d out of %d\n', iAnimals, length(dataPath));
        toc;
    end
end

%% figure1 : Show real data: wiedefield overview - visual trials
cMovie = arrayShrink(nanmean(allData{1},3),mask,'split');
cMovie = cMovie - mean(cMovie(:,:,1:15),3);
cRange = [-0.006 0.006];
% cRange = [-nanmean(abs(cMovie(:))*2) nanmean(abs(cMovie(:)))*2];
figure;
for iSegs = 1:5
    subplot(1,5,iSegs)
    mapImg = imshow(nanmean(cMovie(:,:,segIdxRealign{iSegs+1}),3),[cRange]);
    colormap(mapImg.Parent,colormap_blueblackred); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['All mice. cVision: ' segLabels{iSegs+1}])
    
    hold(mapImg.Parent, 'on');
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end

%% create reconstructed PSTH maps
%non-motor reconstruction
figure
cMovie = (arrayShrink(nanmean(sensoryRec{3},3),mask,'split'));
cRange = [-0.003 0.003];
% cRange = [-nanmean(abs(cMovie(:))*1.5) nanmean(abs(cMovie(:)))*1.5];
for iSegs = 2:4
    subplot(1,3,iSegs-1)
    mapImg = imshow(nanmean(cMovie(:,:,segIdxRealign{iSegs+1}),3),[cRange]);
    colormap(mapImg.Parent,colormap_blueblackred(256)); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['All mice. Sensory rebuild: ' segLabels{iSegs+1}])
    
    hold(mapImg.Parent, 'on');
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end

%% compute modulation in both motor and task regressors
allVar = arrayShrink(nanmean(sum(respVar(:, iMod, ismember(fullRecLabels,'full'), :, iSeg) ,3),4),mask,'split');
motorVar = arrayShrink(nanmean(sum(respVar(:, iMod, ismember(fullRecLabels,'motor'), :, iSeg) ,3),4),mask,'split');
taskVar = arrayShrink(nanmean(sum(respVar(:, iMod, ismember(fullRecLabels,'task'), :, iSeg) ,3),4),mask,'split');

figure
subplot(1,3,1);
mapImg = imshow(taskVar,[0 cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256));hold on; colorbar
title('Mean reconstructed variance - Task');
for x = 1: length(rightHs)
    plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,3,2);
cRange = prctile(motorVar(:),95.9);
mapImg = imshow(motorVar,[0 cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256));hold on; colorbar
title('Mean reconstructed variance - Motor');
for x = 1: length(rightHs)
    plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,3,3); hold off
modIndex = ((taskVar-motorVar) ./ (taskVar + motorVar) + 1) / 2;
mapImg = imshow(modIndex,[0.3 0.7]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256));hold on; colorbar
title('Task index (Motor-Task / Motor+Task)');
hold(mapImg.Parent, 'on');
for x = 1: length(rightHs)
    plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

%% figure 2 : compare real versus fully modeled data
cTitle = {'Data / Model' 'All mice - visual trials'};
mapData = nanmean(nanmean(fullRec{1}(:,segIdxRealign{3},:),2),3);
cRange = [-nanmean(abs(allData{1}(:))*2) nanmean(abs(allData{1}(:)))*2];
Widefield_tracePlot(mapData,areaCoord,{allData{1} fullRec{1}},mask,trialSegments,cRange,colormap_blueblackred,[],cTitle,xRange);

% figure 2: show sensory, motor and cognitive traces
cTitle = {'Full / Motor / Sensory' 'All mice - visual trials'};
mapData = nanmean(nanmean(allData{1}(:,segIdxRealign{3},:),2),3);
cRange = [-nanmean(abs(allData{1}(:))*2) nanmean(abs(allData{1}(:)))*2];
Widefield_tracePlot(mapData,areaCoord,{allData{1} motorRec{3} sensoryRec{3}},mask,trialSegments,cRange,colormap_blueblackred,[],cTitle,xRange,[-0.01 0.01]);

%% figure 2 : compare sensory/motor events in trained/untrained data
cTitle = {'Sensory: Trained / Untrained' 'All mice'};
cData = {cat(3,sensoryRec{1}(:,:,visExp),sensoryRec{2}(:,:,audExp)) cat(3,sensoryRec{2}(:,:,visExp),sensoryRec{1}(:,:,audExp))};
mapData = nanmean(nanmean(cData{1}(:,segIdxRealign{5},:),3),2);
Widefield_tracePlot(mapData,areaCoord,cData,mask,trialSegments,[],colormap_blueblackred,[],cTitle,xRange,[-0.0075 0.01]);

cTitle = {'Motor: Trained / Untrained' 'All mice'};
cData = {cat(3,motorRec{1}(:,:,visExp),motorRec{2}(:,:,audExp)) cat(3,motorRec{2}(:,:,visExp),motorRec{1}(:,:,audExp))};
mapData = nanmean(nanmean(cData{1}(:,segIdxRealign{5},:),3),2);
Widefield_tracePlot(mapData,areaCoord,cData,mask,trialSegments,[],colormap_blueblackred,[],cTitle,xRange,[-0.005 0.02]);

%% figure 2 : compute difference between trained/untrained data for sensory/motor events 
cTitle = {'Sensory: Trained / Untrained' 'All mice'};
cData = {cat(3,sensoryRec{1}(:,:,visExp),sensoryRec{2}(:,:,audExp)) cat(3,sensoryRec{2}(:,:,visExp),sensoryRec{1}(:,:,audExp))};
mapData = nanmean(nanmean(sensoryRec{1}(:,segIdxRealign{4},:),3),2);
Widefield_tracePlot(mapData,areaCoord,cData,mask,trialSegments,[-1 1],colormap_blueblackred,[],cTitle,xRange,[-0.0075 0.0125]);

cTitle = {'Motor: Trained / Untrained' 'All mice'};
cData = {cat(3,motorRec{1}(:,:,visExp),motorRec{2}(:,:,audExp)) cat(3,motorRec{2}(:,:,visExp),motorRec{1}(:,:,audExp))};
mapData = nanmean(nanmean(motorRec{1}(:,segIdxRealign{4},:),3),2);
Widefield_tracePlot(mapData,areaCoord,cData,mask,trialSegments,[-1 1],colormap_blueblackred,[],cTitle,xRange,[-0.0075 0.0125]);

%% figure - show mean amp change and regression coefficients
iMod = 3; %modality for figure. 1 = vision, 2 = audio, 3 = all
modLabel = {'Trained' 'Naive' 'All'}; %labels for modality
maskIdx = find(~mask);

figure;
idx = zeros(1,length(fullRecLabels));
idx(ismember(fullRecLabels,sensorLabels)) = 2;
idx(ismember(fullRecLabels,oMotorLabels)) = 1;
idx(ismember(fullRecLabels,cogLabels)) = 3;

[ax, ampIdx] = regressorPlot(squeeze(nanmean(respVar(:,iMod,idx > 0,:,4)))',fullRecLabels(idx > 0),5,subplot(2,3,4:6),'g',idx(idx > 0),0.6);
title(ax,'PSTH absolute deviation');

%% compute difference between expert reward and novice reward
figure
subplot(1,3,1)
cMap1 = nanmean(squeeze(arrayShrink(respVar(:,2,ismember(fullRecLabels,'visReward'),:,4),mask,'split')),3);
cRange = [0 prctile(cMap1(:),99)];
mapImg = imshow(cMap1,cRange);
colormap(mapImg.Parent,inferno(256));hold on;
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
title('Absolute deviation: Expert reward')

subplot(1,3,2)
cMap2 = nanmean(squeeze(arrayShrink(respVar(:,2,ismember(fullRecLabels,'audReward'),:,4),mask,'split')),3);
cRange = [0 prctile(cMap2(:),95)];
mapImg = imshow(cMap2,cRange);
colormap(mapImg.Parent,inferno(256));hold on;
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
title('Absolute deviation: Novice reward')

% subplot(1,3,3);
% cMap = cMap1 - cMap2;
% cRange = [prctile(cMap(:),10) prctile(cMap(:),90)];
% mapImg = imshow(cMap1 - cMap2,cRange);
% colormap(mapImg.Parent,inferno(256));hold on;
% set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
% title('Expert / Novice reward difference')

hold(mapImg.Parent, 'on');
for x = 1: length(rightHs)
    plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end