%% basic variables
cPath = 'Y:\data\BpodImager\Animals\'; %data path on the server
trialSegments = [0 55 80 130 160]; %segments to create maps for whole trial regressors
motorIdx = 16; %index for zero-lag motor regressor
baseRange = 1:15; %baseline frames
areaIdx = {'VISp' 'RSPd' 'SSp-ll' 'MOs'}; % V1 - RS - HL - M2

%% session data
[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, segIdxRealign] = delayDecRecordings;
sensorLabels = [sensorLabels cogLabels]; %combine sensory and cognitive labels, so sensory is really 'non-motor' or 'task' here
segIdxRealign{2} = 46:65;
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts

%% get reconstructed Vs, modality index, datapath and allOpts
taskV = BpodImager_motorReconstruct(cPath, 'All', 'task', true); taskV = taskV{1}; %get reconstructed V, used full model
opMotorV = BpodImager_motorReconstruct(cPath, 'All', 'opMotor', true); opMotorV = opMotorV{1}; %get reconstructed V, used full model
[spMotorV, recLabels, dataPath, allOpts, modIdx, sideIdx, alignIdx, baseLength, frames, stimTimes] = BpodImager_motorReconstruct(cPath, 'All', 'spMotor', true);  spMotorV = spMotorV{1}; %get reconstructed V, used full model
dimCnt = size(taskV,1);

areaIdx = {'VISp' 'RSPd' 'SSp-ll' 'MOs'}; % V1 - RS - HL - M2
dataPath{1} = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\mSM30\SpatialDisc\10-Oct-2017\';
%get allen maps
load('allenDorsalMapSM.mat')
mask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
[x1, y1] = size(mask);
load([dataPath{1} 'snapshot_1.mat'])
load([dataPath{1} 'opts2.mat'])
snap = alignAllenTransIm(single(snap),opts.transParams);
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

% pre-allocate arrays
allSnap = NaN(sum(~mask(:)), length(dataPath), 'single'); %pre-allocate larger data array for brain picture
allU = NaN(sum(~mask(:)), dimCnt, length(dataPath), 'single'); %pre-allocate larger data array for all Us
allV = NaN(dimCnt, frames, 5, length(dataPath), 'single'); %pre-allocate data array for PSTHs

dSegs = [3 4 5]; % segments from segIdx for d' calculation
dpData{1} = NaN(sum(~mask(:)), length(dSegs), length(dataPath), 'single'); %pre-allocate data array for d' maps
dpData{2} = NaN(sum(~mask(:)), length(dSegs), length(dataPath), 'single'); %pre-allocate data array for d' maps
dpData{3} = NaN(sum(~mask(:)), length(dSegs), length(dataPath), 'single'); %pre-allocate data array for d' maps
dpData{4} = NaN(sum(~mask(:)), length(dSegs), length(dataPath), 'single'); %pre-allocate data array for d' maps

%% load raw data and get aligned U + averaged Vs for vision, audio and all trials
for iAnimals = 1:length(dataPath)

    load([dataPath{iAnimals} 'Vc.mat'],'U')
    load([dataPath{iAnimals} 'snapshot_1.mat'])
    U = U(:,:,1:dimCnt);
    
    snap = alignAllenTransIm(single(snap),allOpts(iAnimals).opts.transParams);
    U = alignAllenTransIm(U,allOpts(iAnimals).opts.transParams);
    
    allSnap(:,iAnimals) = arrayShrink(snap(1:size(mask,1),1:size(mask,2)),mask);
    allU(:,:,iAnimals) = arrayShrink(U(1:size(mask,1),1:size(mask,2),:),mask);
    
    % get Vc that was used for model and get PSTH for correct vis and aud trials. 
    % 3rd and 4th row are correct left and right trials, respectively.
    load([dataPath{iAnimals} 'interpVc.mat'],'Vc')
    Vc = Vc(:,alignIdx{iAnimals});
    Vc = bsxfun(@minus, Vc, mean(Vc,2));
    
    allV(:,:,1,iAnimals) = mean(reshape(Vc(:,modIdx{iAnimals,1}),dimCnt,frames,[]),3); %correct vision
    allV(:,:,2,iAnimals) = mean(reshape(Vc(:,modIdx{iAnimals,2}),dimCnt,frames,[]),3); %correct audio
    allV(:,:,3,iAnimals) = mean(reshape(Vc(:,modIdx{iAnimals,3}),dimCnt,frames,[]),3); %all correct
    allV(:,:,4,iAnimals) = mean(reshape(Vc(:,sideIdx{iAnimals,1} & modIdx{iAnimals,1}),dimCnt,frames,[]),3); %correct vision left
    allV(:,:,5,iAnimals) = mean(reshape(Vc(:,sideIdx{iAnimals,2} & modIdx{iAnimals,1}),dimCnt,frames,[]),3); %correct vision right
    Vc = reshape(Vc, size(Vc,1), frames, []); 
    
    load([dataPath{iAnimals} 'regData.mat'])
    load([dataPath{iAnimals} 'dimBeta.mat'])
    fullR = fullR(alignIdx{iAnimals},:);
    fullR = bsxfun(@minus, fullR, mean(fullR));
    
    cInd = ismember(recIdx(~idx), find(~ismember(recLabels,motorLabels))); %find task
    rVtask = (fullR(:, cInd) * dimBeta(cInd,:))';
    rVtask = reshape(rVtask, size(rVtask,1), frames, []); 

    cIdx = unique(ceil(find(modIdx{iAnimals,3})/frames)); %index for all unisensory, succesful trial
    for iSegs = 1:length(dSegs)
        
        cData = squeeze(mean(Vc(:, segIdxRealign{dSegs(iSegs)},cIdx),2)); %get data for current segment
        covV = cov(cData');  % S x S
        varP = sum((allU(:,:,iAnimals) * covV) .* allU(:,:,iAnimals), 2);  % 1 x P
        stdP = sqrt(varP); %standard deviation map for correct trials
        
        tIdx1 = ismember(cIdx,unique(ceil(find(modIdx{iAnimals,1})/frames))); %index for visual trials
        dpData{1}(:, iSegs, iAnimals) = allU(:,:,iAnimals) * (mean(cData(:,tIdx1),2) - mean(cData(:,~tIdx1),2)); % vision - modality PSTH
        dpData{1}(:, iSegs, iAnimals) = dpData{1}(:, iSegs, iAnimals) ./ stdP; % d' data. modality
        
        tIdx2 = ismember(cIdx,unique(ceil(find(sideIdx{iAnimals,1})/frames))); %index for left trials
        dpData{2}(:, iSegs, iAnimals) = allU(:,:,iAnimals) * (mean(cData(:,tIdx2 & tIdx1),2) - mean(cData(:,~tIdx2 & tIdx1),2)); % left - right PSTH
        dpData{2}(:, iSegs, iAnimals) = dpData{2}(:, iSegs, iAnimals) ./ stdP; % d' data. L/R choice
        
        %determine trained/untrained modality and compute expert vs novice d'
        if visExp(iAnimals)  %visual expert
            dpData{3}(:, iSegs, iAnimals) = allU(:,:,iAnimals) * (mean(cData(:,tIdx1),2) - mean(cData(:,~tIdx1),2)); % vision - modality PSTH
        elseif ~visExp(iAnimals)  %audio expert
            dpData{3}(:, iSegs, iAnimals) = allU(:,:,iAnimals) * (mean(cData(:,~tIdx1),2) - mean(cData(:,tIdx1),2)); % audio - modality PSTH
        end
        dpData{3}(:, iSegs, iAnimals) = dpData{3}(:, iSegs, iAnimals) ./ stdP; % d' data. expertise
        
        %determine trained/untrained modality and compute expert vs novice d' using task reconstruction alone
        cData = squeeze(mean(rVtask(:, segIdxRealign{dSegs(iSegs)},cIdx),2)); %get data for current segment
        covV = cov(cData');  % S x S
        varP = sum((allU(:,:,iAnimals) * covV) .* allU(:,:,iAnimals), 2);  % 1 x P
        stdP = sqrt(varP); %standard deviation map for correct trials
                
        if visExp(iAnimals)  %visual expert
            dpData{4}(:, iSegs, iAnimals) = allU(:,:,iAnimals) * (mean(cData(:,tIdx1),2) - mean(cData(:,~tIdx1),2)); % vision - modality PSTH
        elseif ~visExp(iAnimals)  %audio expert
            dpData{4}(:, iSegs, iAnimals) = allU(:,:,iAnimals) * (mean(cData(:,~tIdx1),2) - mean(cData(:,tIdx1),2)); % audio - modality PSTH
        end
        dpData{4}(:, iSegs, iAnimals) = dpData{4}(:, iSegs, iAnimals) ./ stdP; % d' data. expertise
    
    end
end
clear Vc U dimBeta fullR

%% get reconstructed data
cMod = 1; %use all trials
allData = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
opMotorRec = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
spMotorRec = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
taskRec = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
allTraces = NaN(size(allV,2), length(areaCoord), size(allV,3), size(allV,4));

% convolve real and reconstructed Vs to create PSTHs
for iAnimals = 1:length(dataPath)

    allData(:, :, iAnimals) = (allU(:,:,iAnimals) * allV(:,:,cMod,iAnimals)); % original data
    opMotorRec(:, :, iAnimals) = (allU(:,:,iAnimals) * opMotorV(:,:,cMod,iAnimals)); % operant motor
    spMotorRec(:, :, iAnimals) = (allU(:,:,iAnimals) * spMotorV(:,:,cMod,iAnimals)); % spont. motor
    taskRec(:, :, iAnimals) = (allU(:,:,iAnimals) * taskV(:,:,cMod,iAnimals)); % task
    
    for iAreas = 1 : length(areaCoord)
        traceIdx = ismember(find(~mask),find(areaCoord{iAreas}));
        for iMods = 1 : size(allV,3)
            allTraces(:, iAreas, iMods, iAnimals) = nanmean(allU(traceIdx,:,iAnimals) * allV(:,:,iMods,iAnimals));
        end
    end    
    
end
fullRec = sum(cat(4,taskRec,spMotorRec,opMotorRec),4); %full model reconstruction

%% show dPrime maps
cLabel = {'d" Modality - Stim' 'd" Choice - Stim' 'd" Trained - Stim' 'd" Modality - Delay' 'd" Choice - Delay' 'd" Trained - Delay'};
figure
Cnt = 0;
for iSegs = [1 3]
    for iMod = 1:2
        Cnt = Cnt+1;
        subplot(2,2,Cnt)
        
%         cMap = arrayShrink(nanmean([nanmean(dpData{iMod}(:,iSegs,visExp),3) nanmean(dpData{iMod}(:,iSegs,audExp),3)],2),mask,'split');
        cMap = arrayShrink(nanmean(dpData{iMod}(:,iSegs,:),3),mask,'split');
        
        mapImg = imshow(cMap,[-1 1]);
        colormap(mapImg.Parent,colormap_blueblackred(256)); axis image
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        title(cLabel{Cnt})
        
        hold(mapImg.Parent, 'on');
        for x = 1: length(rightHs)
            plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
        end
    end
end

%%
% cMovie = arrayShrink(nanmean(allData,3),mask,'split');
modLabel = {'Raw' 'Task' 'Spont' 'Instructed'};
cRange = [-0.003 0.003];
    
for iAnimal = 1
% for iAnimal = 1:size(visExp,1)
figure;
% set(gcf,'position',[-1919          41        1920         963]);
Cnt = 0;
% for iMod = 0:3
for iMod = 0
    if iMod == 1
%         cMovie = arrayShrink(nanmean(taskRec(:,:,iAnimal),3),mask,'split');
        cMovie = arrayShrink(nanmean(taskRec,3),mask,'split');
%         cRange = [-nanmean(abs(cMovie(:))*2) nanmean(abs(cMovie(:)))*2];
    elseif iMod == 2
        cMovie = arrayShrink(nanmean(spMotorRec,3),mask,'split');
    elseif iMod == 3
        cMovie = arrayShrink(nanmean(opMotorRec,3),mask,'split');        
    elseif iMod == 0
        cMovie = arrayShrink(nanmean(allData,3),mask,'split');        
    end
    cMovie = cMovie - mean(cMovie(:,:,1:15),3);
    if iMod == 0
    cRange = [-nanmean(abs(cMovie(:))*2.5) nanmean(abs(cMovie(:)))*2.5];
    end

    for iSegs = 1:5
        Cnt = Cnt + 1;
        subplot(4,5,Cnt)
        mapImg = imshow(nanmean(cMovie(:,:,segIdxRealign{iSegs+1}),3),[cRange]);
        colormap(mapImg.Parent,colormap_blueblackred(256)); axis image
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        title([modLabel{iMod+1} ' - ' segLabels{iSegs+1}])
        
        hold(mapImg.Parent, 'on');
        for x = 1: length(rightHs)
            plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
        end
        drawnow;

    end
end
% Widefield_tracePlot(nanmean(taskRec(:,85,iAnimal),3),areaCoord,{allData(:,:,iAnimal), spMotorRec(:,:,iAnimal)},mask,trialSegments,cRange,colormap_blueblackred,[],['iAnimal: ' num2str(iAnimal)],xRange,[cRange(1)*2 cRange(2)*2], false, true);
% Widefield_tracePlot(nanmean(taskRec(:,85,:),3),areaCoord,{allData},mask,trialSegments,cRange,colormap_blueblackred,[],['iAnimal: ' num2str(iAnimal)],xRange,[cRange(1)*3 cRange(2)*4.5]);
drawnow;
end

%% make traces figures
baseLength = 79;
% vision - audio figure
figure; 

for iAreas = 1 : length(areaCoord)
   
    subplot(length(areaCoord), 1, iAreas);
    cData = squeeze(allTraces(:,iAreas,1,:))';
    cData = cData - mean(cData(:,1:15),2);
    timeAx = 1/30 : 1/30 : size(cData,2)/30;
    stdshade(cData(:, 1:baseLength), 'k', timeAx(1:baseLength), 0.5); hold on
    stdshade(cData(:, baseLength+1:end), 'k', timeAx(baseLength+1:end), 0.5);
    
    cData = squeeze(allTraces(:,iAreas,2,:))';
    cData = cData - mean(cData(:,1:15),2);
    timeAx = 1/30 : 1/30 : size(cData,2)/30;
    stdshade(cData(:, 1:baseLength), 'r', timeAx(1:baseLength), 0.5); hold on
    stdshade(cData(:, baseLength+1:end), 'r', timeAx(baseLength+1:end), 0.5);
    
    xlim([1/30 size(allTraces,1) / 30]); ylim([-0.01 0.02]);
    title([areaIdx{iAreas} ' - Vision / audio responses']);
    vline(([55 baseLength + [1 19 34 52 82]]) ./ 30); hline(0);
    
end

%% vision - left/right vision figure
figure; 

for iAreas = 1 : length(areaCoord)
   
    subplot(length(areaCoord), 1, iAreas);
    cData = squeeze(allTraces(:,iAreas,4,:))';
    cData = cData - mean(cData(:,1:15),2);
    timeAx = 1/30 : 1/30 : size(cData,2)/30;
    stdshade(cData(:, 1:baseLength), 'k', timeAx(1:baseLength), 0.5); hold on
    stdshade(cData(:, baseLength+1:end), 'k', timeAx(baseLength+1:end), 0.5);
    
    cData = squeeze(allTraces(:,iAreas,5,:))';
    cData = cData - mean(cData(:,1:15),2);
    timeAx = 1/30 : 1/30 : size(cData,2)/30;
    stdshade(cData(:, 1:baseLength), 'r', timeAx(1:baseLength), 0.5); hold on
    stdshade(cData(:, baseLength+1:end), 'r', timeAx(baseLength+1:end), 0.5);
    
    xlim([1/30 size(allTraces,1) / 30]); ylim([-0.01 0.02]);
    title([areaIdx{iAreas} ' - Left / right vision responses']);
    vline(([55 baseLength + [1 19 34 52 82]]) ./ 30); hline(0);
    
end

%% exp. vs modality figure
figure; 

for iAreas = 1 : length(areaCoord)
   
    subplot(length(areaCoord), 1, iAreas);
    cData = squeeze(cat(4,allTraces(:,iAreas,1,visExp),allTraces(:,iAreas,2,audExp)))';
    cData = cData - mean(cData(:,1:15),2);
    timeAx = 1/30 : 1/30 : size(cData,2)/30;
    stdshade(cData(:, 1:baseLength), 'k', timeAx(1:baseLength), 0.5); hold on
    stdshade(cData(:, baseLength+1:end), 'k', timeAx(baseLength+1:end), 0.5);
    
    cData = squeeze(cat(4,allTraces(:,iAreas,2,visExp),allTraces(:,iAreas,1,audExp)))';
    cData = cData - mean(cData(:,1:15),2);
    timeAx = 1/30 : 1/30 : size(cData,2)/30;
    stdshade(cData(:, 1:baseLength), 'r', timeAx(1:baseLength), 0.5); hold on
    stdshade(cData(:, baseLength+1:end), 'r', timeAx(baseLength+1:end), 0.5);
    
    xlim([1/30 size(allTraces,1) / 30]); 
%     ylim([-0.01 0.02]);
    title([areaIdx{iAreas} ' - Exp / novice modality responses']);
    vline(([55 baseLength + [1 19 34 52 82]]) ./ 30); hline(0);
    
end