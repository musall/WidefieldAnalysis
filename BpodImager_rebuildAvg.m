%% basic variables
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
trialSegments = [0 55 80 130 160]; %segments to create maps for whole trial regressors
motorIdx = 16; %index for zero-lag motor regressor
baseRange = 1:15; %baseline frames
areaIdx = {'VISp' 'RSPd' 'SSp-ll' 'MOs'}; % V1 - RS - HL - M2
normData = false; %flag to zscore data

%% session data
[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, segIdxRealign] = delayDecRecordings;
sensorLabels = [sensorLabels cogLabels]; %combine sensory and cognitive labels, so sensory is really 'non-motor' or 'task' here
segIdxRealign{2} = 46:65;
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts

%% get reconstructed Vs, modality index, datapath and allOpts
[recV, recLabels, dataPath, allOpts, modIdx, sideIdx, alignIdx, baseLength, frames, stimTimes] = BpodImager_motorReconstruct(cPath, 'All', 'All', false); %get reconstructed V, used full model
recV = cat(5,recV{:});
dimCnt = size(recV,1);
stimTimes = cat(2,stimTimes{:});
stimTimes(stimTimes > 1.5) = [];
figure; histogram(stimTimes,50); axis square; %show stim time histogram
nRecLabels = [{'full'} recLabels{:} {'motor' 'motorVideo' 'sensory' 'cognitive'}];

 %get allen maps
load('allenDorsalMapSM.mat')
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
allSnap = NaN(sum(~mask(:)), length(dataPath), 'single'); %pre-allocate larger data array for brain picture
allU = NaN(sum(~mask(:)), dimCnt, length(dataPath), 'single'); %pre-allocate larger data array for all Us
allV = NaN(dimCnt, frames, 4, length(dataPath), 'single'); %pre-allocate data array for PSTHs
if normData
    allStd = NaN(sum(~mask(:)), length(dataPath), 'single'); %pre-allocate larger data array for brain picture
    modelStd = NaN(sum(~mask(:)), length(dataPath), 'single'); %pre-allocate larger data array for brain picture
else
    allStd = ones(1,length(dataPath)); %std maps are not needed if no zscoring is done
    modelStd = ones(1,length(dataPath)); %std maps are not needed if no zscoring is done
end

%% load raw data
for iAnimals = 1:length(dataPath)

    load([dataPath{iAnimals} 'Vc.mat'],'U')
    load([dataPath{iAnimals} 'snapshot_1.mat'])
    U = U(:,:,1:dimCnt);
    
    snap = alignAllenTransIm(single(snap),allOpts(iAnimals).opts.transParams);
    U = alignAllenTransIm(U,allOpts(iAnimals).opts.transParams);
    
    allSnap(:,iAnimals) = arrayShrink(snap(1:size(mask,1),1:size(mask,2)),mask);
    allU(:,:,iAnimals) = arrayShrink(U(1:size(mask,1),1:size(mask,2),:),mask);
    
    % get Vc that was used for model and get PSTH for correct vis and aud
    % trials. 3rd and 4th row are correct left and right trials, respectively.
    load([dataPath{iAnimals} 'interpVc.mat'],'Vc')
    Vc = Vc(:,alignIdx{iAnimals});
    
    allV(:,:,1,iAnimals) = mean(reshape(Vc(:,modIdx{iAnimals,1}),dimCnt,frames,[]),3);
    allV(:,:,2,iAnimals) = mean(reshape(Vc(:,modIdx{iAnimals,2}),dimCnt,frames,[]),3);

    allV(:,:,3,iAnimals) = mean(reshape(Vc(:,sideIdx{iAnimals,1} & modIdx{iAnimals,1}),dimCnt,frames,[]),3); %left vision trials
    allV(:,:,4,iAnimals) = mean(reshape(Vc(:,sideIdx{iAnimals,2} & modIdx{iAnimals,1}),dimCnt,frames,[]),3); %right vision trials

    if normData
        % get standard deviation over all succesful trials
        covV = cov(Vc(:,modIdx{iAnimals,3})');  % S x S
        allStd(:, iAnimals) = sqrt(sum((squeeze(allU(:,:,iAnimals)) * covV) .* squeeze(allU(:,:,iAnimals)), 2));  % 1 x P

        % get standard deviation over all modeled succesful trials
        load([dataPath{iAnimals} 'dimBeta.mat'],'dimBeta')
        load([dataPath{iAnimals} 'regData.mat'],'fullR'); %load model data
        
        covV = cov(fullR(modIdx{iAnimals,3},:) * dimBeta);  % covariance of full model
        modelStd(:, iAnimals) = sqrt(sum((squeeze(allU(:,:,iAnimals)) * covV) .* squeeze(allU(:,:,iAnimals)), 2));  % 1 x P
        clear Vc U dimBeta fullR
    end
end

%% get reconstructed data
fullRec{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
fullRec{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

motorRec{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
motorRec{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

sensoryRec{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
sensoryRec{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

allData{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
allData{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

sideData{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
sideData{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

dpData{1} = NaN(sum(~mask(:)), 2, length(dataPath), 'single'); %pre-allocate data array for PSTHs
dpData{2} = NaN(sum(~mask(:)), 2, length(dataPath), 'single'); %pre-allocate data array for PSTHs
dpData{3} = NaN(sum(~mask(:)), 2, length(dataPath), 'single'); %pre-allocate data array for PSTHs

% frame index for stimulus and waiting segments
cSegs = cat(2,segIdxRealign{ismember(segLabels,{'Stim1' 'Stim2'})}); %time index for stimulus segments
cSegs = [cSegs(1) segIdxRealign{ismember(segLabels,{'Wait'})}(round(end/2))]; %add delay period
cSegs = cSegs([1 end]);
respAngle = NaN(sum(~mask(:)), 3, size(recV,5) + 3, length(dataPath), 'single');
respVar = NaN(sum(~mask(:)), 3, size(recV,5) + 3, length(dataPath), 'single');

% convolve real and reconstructed Vs to create PSTHs
for iAnimals = 1:length(dataPath)

    allData{1}(:, :, iAnimals) = (allU(:,:,iAnimals) * allV(:,:,1,iAnimals)) ./ allStd(:, iAnimals); % original data
    allData{2}(:, :, iAnimals) = (allU(:,:,iAnimals) * allV(:,:,2,iAnimals)) ./ allStd(:, iAnimals); % original data

    sideData{1}(:, :, iAnimals) = allU(:,:,iAnimals) * allV(:,:,3,iAnimals) ./ allStd(:, iAnimals); % original data
    sideData{2}(:, :, iAnimals) = allU(:,:,iAnimals) * allV(:,:,4,iAnimals) ./ allStd(:, iAnimals); % original data
         
    fullRec{1}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,1,iAnimals,:),5) ./ allStd(:, iAnimals); % full reconstruction
    fullRec{2}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,2,iAnimals,:),5) ./ allStd(:, iAnimals); % full reconstruction
    
    motorRec{1}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,1,iAnimals,ismember(recLabels,motorLabels)),5) ./ allStd(:, iAnimals); % motor reconstruction
    motorRec{2}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,2,iAnimals,ismember(recLabels,motorLabels)),5) ./ allStd(:, iAnimals); % motor reconstruction
    
    sensoryRec{1}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,1,iAnimals,ismember(recLabels,sensorLabels)),5) ./ allStd(:, iAnimals); % non-motor reconstruction
    sensoryRec{2}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,2,iAnimals,ismember(recLabels,sensorLabels)),5) ./ allStd(:, iAnimals); % non-motor reconstruction
    
    % reload Vc to compute standard deviation later
    load([dataPath{iAnimals} 'interpVc.mat'],'Vc')
    Vc = Vc(:,alignIdx{iAnimals});
    Vc = reshape(Vc,size(Vc,1),frames,[]);
    
    x = [3 5]; % segments from segIdx for d' calculation
    cIdx = unique(ceil(find(modIdx{iAnimals,3})/frames)); %index for all unisensory, succesful trial
    
    for iSegs = 1:2
        % compare standard deviation in current trial segment
        cData = squeeze(mean(Vc(:, segIdxRealign{x(iSegs)},cIdx),2));
        covV = cov(cData');  % S x S
        varP = sum((allU(:,:,iAnimals) * covV) .* allU(:,:,iAnimals), 2);  % 1 x P
        stdP = sqrt(varP); %standard deviation map for correct trials
        
        tIdx1 = ismember(cIdx,unique(ceil(find(modIdx{iAnimals,1})/frames))); %index for visual trials
        dpData{1}(:, iSegs, iAnimals) = allU(:,:,iAnimals) * (mean(cData(:,tIdx1),2) - mean(cData(:,~tIdx1),2)); % vision - modality PSTH
        dpData{1}(:, iSegs, iAnimals) = dpData{1}(:, iSegs, iAnimals) ./ stdP; % d' data. modality
        
        tIdx2 = ismember(cIdx,unique(ceil(find(sideIdx{iAnimals,1})/frames))); %index for left trials
        dpData{2}(:, iSegs, iAnimals) = allU(:,:,iAnimals) * (mean(cData(:,tIdx2 & tIdx1),2) - mean(cData(:,~tIdx2 & tIdx1),2)); % left - right PSTH
        dpData{2}(:, iSegs, iAnimals) = dpData{2}(:, iSegs, iAnimals) ./ stdP; % d' data. L/R choice
        
        %determine trained/untrained modality for iMods 1/2.
        if visExp(iAnimals)  %visual expert
            dpData{3}(:, iSegs, iAnimals) = allU(:,:,iAnimals) * (mean(cData(:,tIdx1),2) - mean(cData(:,~tIdx1),2)); % vision - modality PSTH
        elseif ~visExp(iAnimals)  %audio expert
            dpData{3}(:, iSegs, iAnimals) = allU(:,:,iAnimals) * (mean(cData(:,~tIdx1),2) - mean(cData(:,tIdx1),2)); % audio - modality PSTH
        end
        dpData{3}(:, iSegs, iAnimals) = dpData{3}(:, iSegs, iAnimals) ./ stdP; % d' data. expertise
    end
    
    % Check for ramp angle during stimulus and waiting period in each regressor
    for iMods = 1:3
%     for iMods = 3
        for iRegs = 1 : size(recV,5) + 3
%         for iRegs = size(recV,5) + 1 : size(recV,5) + 3
            if iRegs <= size(recV,5)
                cIdx = iRegs;
            elseif iRegs == size(recV,5) + 1
                cIdx = true(1,length(recLabels)); %all regressors
            elseif iRegs == size(recV,5) + 2
                cIdx = ismember(recLabels, motorLabels); %motor regressors
            elseif iRegs == size(recV,5) + 3
                cIdx = ~ismember(recLabels, motorLabels); %task regressors
            end
            
             %determine trained/untrained modality for iMods 1/2. 3 is full PSTH.
            if iMods == 1 && visExp(iAnimals) || iMods == 2 && ~visExp(iAnimals)
                cMovie = allU(:,:,iAnimals) * sum(recV(:,:,1,iAnimals,cIdx),5);
            elseif iMods == 1 && ~visExp(iAnimals) || iMods == 2 && visExp(iAnimals)
                cMovie = allU(:,:,iAnimals) * sum(recV(:,:,2,iAnimals,cIdx),5);
            elseif iMods == 3
                cMovie = allU(:,:,iAnimals) * sum(recV(:,:,3,iAnimals,cIdx),5);
            end
            
            cBetas = [ones(diff(cSegs)+1,1) (1 : diff(cSegs)+1)'] \ (cMovie(:,cSegs(1):cSegs(2))');
            respAngle(:,iMods,iRegs,iAnimals) = cBetas(2,:)'; %slope during stimulus period
            respVar(:,iMods,iRegs,iAnimals) = mean(abs(cMovie(:,cSegs(1):cSegs(2))),2); %signal deviation from zero during stimulus/delay period
%             respVar(:,iMods,iRegs,iAnimals) = sum(abs(cMovie(:,:)),2); %signal deviation from zero during stimulus/delay period
            clear cMovie cBetas
        end
    end
end

%% figure1 : Show real data: wiedefield overview - visual trials
cMovie = arrayShrink(nanmean(allData{1},3),mask,'split');
cMovie = cMovie - mean(cMovie(:,:,1:15),3);
cRange = [-nanmean(abs(cMovie(:))*2) nanmean(abs(cMovie(:)))*2];
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

%% figure1 : show real data: vision and auditory
% compare modalities - visual experts mice
cTitle = {'Vision / Audio' 'Visual experts'};
mapData = nanmean(nanmean(allData{1}(:,segIdxRealign{3},:),2),3);
Widefield_tracePlot(mapData,areaCoord,allData,mask,trialSegments,cRange,colormap_blueblackred,[],cTitle,xRange,[cRange(1) cRange(2)*3.5], true);
    
% figure1 : show real data: left and right
cTitle = {'Left / Right' 'Visual experts'};
mapData = nanmean(nanmean(sideData{1}(:,segIdxRealign{3},:),2),3);
Widefield_tracePlot(mapData,areaCoord,sideData,mask,trialSegments,cRange,colormap_blueblackred,[],cTitle,xRange,[cRange(1) cRange(2)*3.5], true);

% figure 1: show real data: trained vs untrained
cTitle = {'Trained / Untrained' 'All mice'};
cData = {cat(3,allData{1}(:,:,visExp),allData{2}(:,:,audExp)) cat(3,allData{2}(:,:,visExp),allData{1}(:,:,audExp))};
mapData = nanmean(nanmean(cData{1}(:,segIdxRealign{3},:),3),2);
Widefield_tracePlot(mapData,areaCoord,cData,mask,trialSegments,cRange,colormap_blueblackred,[],cTitle,xRange,[cRange(1) cRange(2)*3.5], true);


%% figure1 : show d' maps: modality and choice
cLabel = {'d" Modality - Stim' 'd" Choice - Stim' 'd" Trained - Stim' 'd" Modality - Delay' 'd" Choice - Delay' 'd" Trained - Delay'};
figure
Cnt = 0;
for iSegs = 1:2
    for iMod = 1:3
        Cnt = Cnt+1;
        subplot(2,3,Cnt)
        
        cMap = arrayShrink(nanmean(dpData{iMod}(:,iSegs,:),3),mask,'split');
        
        mapImg = imshow(cMap,[-0.5 0.5]);
        colormap(mapImg.Parent,colormap_blueblackred); axis image
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        title(cLabel{Cnt})
        
        hold(mapImg.Parent, 'on');
        for x = 1: length(dorsalMaps.edgeOutlineSplit)
            plot(dorsalMaps.edgeOutlineSplit{x}(:,2), dorsalMaps.edgeOutlineSplit{(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
%             plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
        end
    end
end

%% create reconstructed PSTH maps
%non-motor reconstruction
figure
cMovie = (arrayShrink(nanmean(motorRec{1},3),mask,'split'));
% cRange = [-nanmean(abs(cMovie(:))*2.5) nanmean(abs(cMovie(:)))*2.5];
for iSegs = 1:5
    subplot(1,5,iSegs)
    mapImg = imshow(nanmean(cMovie(:,:,segIdxRealign{iSegs+1}),3),[cRange]);
    colormap(mapImg.Parent,colormap_blueblackred); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['All mice. Sensory rebuild: ' segLabels{iSegs+1}])
    
    hold(mapImg.Parent, 'on');
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end

%% compute modulation in both motor and task regressors
iMod = 3; % 3 = all trials

idx = zeros(1,length(nRecLabels));
idx(ismember(nRecLabels,sensorLabels)) = 1;
idx(ismember(nRecLabels,motorLabels)) = 2;

allVar = arrayShrink(nanmean(sum(respVar(:,iMod,length(recLabels)+1,:) ,3),4),mask,'split');
motorVar = arrayShrink(nanmean(sum(respVar(:,iMod,length(recLabels)+2,:) ,3),4),mask,'split');
taskVar = arrayShrink(nanmean(sum(respVar(:,iMod,length(recLabels)+3,:) ,3),4),mask,'split');
cRange = prctile(allVar(:),95);

figure
subplot(1,4,1);
mapImg = imshow(allVar,[0 cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,jet(256));hold on; colorbar
title('Mean reconstructed variance - All');
for x = 1: length(rightHs)
    plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,4,2);
mapImg = imshow(motorVar,[0 cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,jet(256));hold on; colorbar
title('Mean reconstructed variance - Motor');
for x = 1: length(rightHs)
    plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,4,3);
mapImg = imshow(taskVar,[0 cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,jet(256));hold on; colorbar
title('Mean reconstructed variance - Task');
for x = 1: length(rightHs)
    plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,4,4); hold off
modIndex = ((taskVar-motorVar) ./ (taskVar + motorVar) + 1) / 2;
cRange = prctile(modIndex(:),99);
mapImg = imshow(modIndex,[0.2 0.8]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,jet(256));hold on; colorbar
title('Task index (Motor-Task / Motor+Task)');
hold(mapImg.Parent, 'on');
for x = 1: length(rightHs)
    plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end


%% modality and variance modulation indices
figure

subplot(1,3,2); hold off
cRange = prctile(allVar(:),99);
mapImg = imshow(allVar,[0 cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,'jet');hold on; colorbar
title('All variance');
hold(mapImg.Parent, 'on');
for x = 1: length(rightHs)
    plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,3,3); hold off
varIndex = modIndex .* allVar;
% varIndex =(1-abs(modIndex)) .* allVar;
cRange = prctile(varIndex(:),99);
mapImg = imshow(varIndex,[0 cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,'jet');hold on; colorbar
title('Adjusted task index (TI * Variance)');
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
Widefield_tracePlot(mapData,areaCoord,{fullRec{1} motorRec{1} sensoryRec{1}},mask,trialSegments,cRange,colormap_blueblackred,[],cTitle,xRange,[]);

% figure 2: show sensory traces
cTitle = {'Data / Motor / Sensory' 'All mice - visual trials'};
mapData = nanmean(nanmean(allData{1}(:,segIdxRealign{3},:),2),3);
Widefield_tracePlot(mapData,areaCoord,sensoryRec,mask,trialSegments,cRange,colormap_blueblackred,visExp,cTitle,xRange,[]);

%% figure 2 : model similiarity index
cMovie = allData{1};
cMovie1 = arrayShrink(nanmean(abs(motorRec{1}-cMovie),3),mask,'split');
cMovie2 = arrayShrink(nanmean(abs(sensoryRec{1}-cMovie),3),mask,'split');
cMovie3 = (cMovie2 - cMovie1) ./ (cMovie1 + cMovie2);

figure;
for iSegs = 1:6
    subplot(2,3,iSegs)
    if iSegs ~= 6
        mapImg = imshow(nanmean(cMovie3(:,:,segIdxRealign{iSegs+1}),3),[-0.5 0.5]);
        title(['All mice. Motor rebuild: ' segLabels{iSegs+1}])
    else
        mapImg = imshow(nanmean(cMovie3,3),[-0.5 0.5]);
        title('All mice. Motor rebuild: All')
    end
    colormap(mapImg.Parent,'parula'); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    
    hold(mapImg.Parent, 'on');
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end

%% figure 2 : model component correlation
cMap1 = NaN(size(allData{1},1),1);
cMap2 = NaN(size(allData{1},1),1);
a = nanmean(allData{1},3);
b = nanmean(motorRec{1},3);
c = nanmean(sensoryRec{1},3);

for x = 1:size(allData{1},1)
    cMap1(x) = corr2(a(x,:),b(x,:));
    cMap2(x) = corr2(a(x,:),c(x,:));
end

figure
subplot(1,3,1);
imagesc(arrayShrink(cMap1,mask,'split'));axis square;caxis([0 1]);colorbar;colormap jet
subplot(1,3,2);
imagesc(arrayShrink(cMap2,mask,'split'));axis square;caxis([0 1]);colorbar;colormap jet
subplot(1,3,3); 
cImg = imagesc(arrayShrink(cMap1,mask,'split'));axis square;caxis([0 1])
set(cImg,'AlphaData',arrayShrink(cMap1,mask,'split'));hold on;
cImg = imagesc(arrayShrink(cMap2+1,mask,'split'));axis square;caxis([0 1])
set(cImg,'AlphaData',arrayShrink(cMap1,mask,'split'));
caxis([0 2]);
colormap(cImg.Parent,[[(1:128)/128 zeros(1,128)]; zeros(1,256);[zeros(1,128) ((129:256)/128)-1]]');


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
cData = {cat(3,sensoryRec{1}(:,:,visExp)-sensoryRec{2}(:,:,visExp),sensoryRec{2}(:,:,audExp)-sensoryRec{1}(:,:,audExp))};
mapData = nanmean(nanmean(sensoryRec{1}(:,segIdxRealign{4},:),3),2);
Widefield_tracePlot(mapData,areaCoord,cData,mask,trialSegments,[-1 1],colormap_blueblackred,[],cTitle,xRange,[-0.0075 0.0125]);

cTitle = {'Motor: Trained / Untrained' 'All mice'};
cData = {cat(3,motorRec{1}(:,:,visExp)-motorRec{2}(:,:,visExp),motorRec{2}(:,:,audExp)-motorRec{1}(:,:,audExp))};
mapData = nanmean(nanmean(motorRec{1}(:,segIdxRealign{4},:),3),2);
Widefield_tracePlot(mapData,areaCoord,cData,mask,trialSegments,[-1 1],colormap_blueblackred,[],cTitle,xRange,[-0.0075 0.0125]);


%% compare real data, versus motor and 'sensory' component
modality = {'Vision' 'audio'}; %labels for modality
for iMod = 1:2
    %visual trials
    cTitle = {'Motor / Sensory / Time' ['All mice - ' modality{iMod} ' trials']};
    mapData = nanmean(nanmean(allData{iMod}(:,segIdxRealign{3},:),3),2);
    Widefield_tracePlot(mapData,areaCoord,{motorRec{iMod} sensoryRec{iMod}},mask,trialSegments,[-2 2],colormap_blueblackred,[],cTitle,xRange);
    
end

%% compare motor components across modalities
cTitle = {'Motor: Vision / Audio' 'All mice'};
Widefield_tracePlot(mapData,areaCoord,motorRec,mask,trialSegments,[-2 2],colormap_blueblackred,visExp,cTitle,xRange);

% cTitle = {'Motor: Vision / Audio' 'Visual experts'};
% mapData = nanmean(nanmean(motorRec{1}(:,segIdxRealign{3},visExp),3),2);
% Widefield_tracePlot(mapData,areaCoord,motorRec,mask,trialSegments,[-0.005 0.005],colormap_blueblackred,visExp,cTitle,xRange);
% 
% cTitle = {'Motor: Vision / Audio' 'Auditory experts'};
% mapData = nanmean(nanmean(motorRec{1}(:,segIdxRealign{3},audExp),3),2);
% Widefield_tracePlot(mapData,areaCoord,motorRec,mask,trialSegments,[-0.005 0.005],colormap_blueblackred,audExp,cTitle,xRange);


%% compare sensory components across modalities
cTitle = {'Sensory: Vision / Audio' 'All mice'};
Widefield_tracePlot(mapData,areaCoord,sensoryRec,mask,trialSegments,[-2 2],colormap_blueblackred,visExp,cTitle,xRange);

% cTitle = {'Sensory: Vision / Audio' 'Visual experts'};
% mapData = nanmean(nanmean(sensoryRec{1}(:,segIdxRealign{3},visExp),3),2);
% Widefield_tracePlot(mapData,areaCoord,sensoryRec,mask,trialSegments,[],colormap_blueblackred,visExp,cTitle,xRange);
% 
% cTitle = {'Sensory: Vision / Audio' 'Auditory experts'};
% mapData = nanmean(nanmean(sensoryRec{1}(:,segIdxRealign{3},audExp),3),2);
% Widefield_tracePlot(mapData,areaCoord,sensoryRec,mask,trialSegments,[],colormap_blueblackred,audExp,cTitle,xRange);

%% compute single regressor contribution
iReg = 1; %ID for regressor that should be reconstructed

cRec{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
cRec{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

for iAnimals = 1:length(dataPath)
    cRec{1}(:, :, iAnimals) = allU(:,:,iAnimals) * recV(:,:,1,iAnimals,iReg); % reconstruction
    cRec{2}(:, :, iAnimals) = allU(:,:,iAnimals) * recV(:,:,2,iAnimals,iReg); % reconstruction
end

cTitle = {['Data / Model / ' recLabels{iReg}] 'All mice - visual trials'};
Widefield_tracePlot(mapData,areaCoord,{allData{1} fullRec{1} cRec{1}},mask,trialSegments,[-0.005 0.005],colormap_blueblackred,[],cTitle,xRange);


%% compare trained vs untrained modality. Motor and sensory regressors.
cTitle = {'Sensory: Trained / Untrained' 'All mice'};
cData = {zscore(cat(3,sensoryRec{1}(:,:,visExp),sensoryRec{2}(:,:,audExp)),[],2) zscore(cat(3,sensoryRec{2}(:,:,visExp),sensoryRec{1}(:,:,audExp)),[],2)};
mapData = nanmean(nanmean(cData{1}(:,segIdxRealign{3},:),3),2);
Widefield_tracePlot(mapData,areaCoord,cData,mask,trialSegments,[-1 1],colormap_blueblackred,[],cTitle,xRange,[-2 4]);

% compute modulation index for sensory component
sensoryMod = NaN(size(cData{1},1),length(segIdxRealign),size(cData{1},3));
sensoryDiff = NaN(size(cData{1},1),length(segIdxRealign),size(cData{1},3));
h1 = figure; h2 = figure;
for iSegs = 1 : length(segIdxRealign)
    sensoryMod(:,iSegs,:) = nanmean(cData{1}(:,segIdxRealign{iSegs},:),2);
    sensoryDiff(:,iSegs,:) = nanmean(cData{1}(:,segIdxRealign{iSegs},:),2) - nanmean(cData{2}(:,segIdxRealign{iSegs},:),2);

    figure(h1);
    subplot(2,3,iSegs);
    mapImg = imshow(squeeze(arrayShrink(nanmean(sensoryDiff(:,iSegs,:),3),mask,'split')),[-.75 .75]);
    colormap(mapImg.Parent,colormap_blueblackred); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['Sensory difference: ' segLabels{iSegs}])

    figure(h2);
    subplot(2,3,iSegs);
    mapImg = imshow(squeeze(arrayShrink(nanmean(sensoryMod(:,iSegs,:),3),mask,'split')),[-.75 .75]);
    colormap(mapImg.Parent,colormap_blueblackred); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['Sensory (trained): ' segLabels{iSegs}])
end

cTitle = {'Motor: Trained / Untrained' 'All mice'};
cData = {zscore(cat(3,motorRec{1}(:,:,visExp),motorRec{2}(:,:,audExp)),[],2) zscore(cat(3,motorRec{2}(:,:,visExp),motorRec{1}(:,:,audExp)),[],2)};
mapData = nanmean(nanmean(cData{1}(:,segIdxRealign{3},:),3),2);
Widefield_tracePlot(mapData,areaCoord,cData,mask,trialSegments,[-1 1],colormap_blueblackred,[],cTitle,xRange,[-2 4]);

% compute modulation index for motor component
motorMod = NaN(size(cData{1},1),length(segIdxRealign),size(cData{1},3));
motorDiff = NaN(size(cData{1},1),length(segIdxRealign),size(cData{1},3));
h1 = figure; h2 = figure;
for iSegs = 1 : length(segIdxRealign)
    motorMod(:,iSegs,:) = nanmean(cData{1}(:,segIdxRealign{iSegs},:),2);
    motorDiff(:,iSegs,:) = nanmean(cData{1}(:,segIdxRealign{iSegs},:),2) - nanmean(cData{2}(:,segIdxRealign{iSegs},:),2);
    
    figure(h1);
    subplot(2,3,iSegs);
    mapImg = imshow(squeeze(arrayShrink(nanmean(motorDiff(:,iSegs,:),3),mask,'split')),[-.75 .75]);
    colormap(mapImg.Parent,colormap_blueblackred); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['Motor difference: ' segLabels{iSegs}])

    figure(h2);
    subplot(2,3,iSegs);
    mapImg = imshow(squeeze(arrayShrink(nanmean(motorMod(:,iSegs,:),3),mask,'split')),[-.75 .75]);
    colormap(mapImg.Parent,colormap_blueblackred); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['Motor (trained): ' segLabels{iSegs}])
end
clear cData


%% figure - show mean amp change and regression coefficients
iMod = 1; %modality for figure. 1 = vision, 2 = audio, 3 = all
modLabel = {'Trained' 'Naive' 'All'}; %labels for modality
maskIdx = find(~mask);

% mean angles
figure;
subplot(2,3,1);
cMap = arrayShrink(nanmean(respAngle(:,iMod,size(recV,5) + 1,:),4),mask,'split');
cRange = [-nanmedian(abs(cMap(:))*2) nanmedian(abs(cMap(:)))*2];
mapImg = imshow(cMap,cRange);
set(mapImg,'AlphaData',~mask); %make NaNs transparent.
colormap(mapImg.Parent,colormap_blueblackred);hold on; colorbar
title(['Mean angle - ' modLabel{iMod}]);
traceColor = get(gca,'ColorOrder');
hold(mapImg.Parent, 'on');
for x = ((length(dorsalMaps.edgeOutlineSplit)-1) / 2 + 1):length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{x}(1,:), dorsalMaps.edgeOutlineSplit{x}(2,:), 'w', 'LineWidth', 1);axis image
end
    
cInd = ~ismember(nRecLabels,{'motor' 'motorVideo' 'sensory' 'cognitive'});
idx = zeros(1,length(nRecLabels));
idx(ismember(nRecLabels,sensorLabels)) = 1;
idx(ismember(nRecLabels,motorLabels)) = 2;
idx(ismember(nRecLabels,cogLabels)) = 3;
regressorPlot(squeeze(nanmean(respAngle(:,iMod,:,:)))',{recLabels{:} 'full' 'motor' 'sensory' 'time'},5,subplot(2,3,2:3));


[ax, angleIdx] = regressorPlot(squeeze(nanmean(respAngle(:,iMod,:,:)))',{recLabels{:} 'full' 'motor' 'sensory' 'time'},5,subplot(2,3,2:3));
title(ax,['PSTH angle during stimulus/wait period (' modLabel{iMod} ' trials)']);
hline(0,'--k');

for x = 1 : size(areaCoord,2)
    
    areas(x) = rectangle('Position',[areaCoord(1,x)-areaCoord(3,x) areaCoord(2,x)-areaCoord(3,x) areaCoord(3,x)*2 areaCoord(3,x)*2],'Curvature',1,'linewidth',2,'linestyle','--','Edgecolor','w','parent',mapImg.Parent);
    text(areaCoord(1,x)-6,areaCoord(2,x),num2str(x),'color','w','FontSize',15,'parent',mapImg.Parent)
    
    traceIdx{x} = find(createCirclesMask(mapImg.CData,(areas(x).Position(1:2) + (areas(x).Position(3)/2)),areas(x).Position(3)/2)); %get mask from current area and find index
    traceIdx{x} = ismember(maskIdx,traceIdx{x});
    
    ax = regressorPlot(squeeze(nanmean(respAngle(traceIdx{x},iMod,:,:)))',{recLabels{:} 'full' 'motor' 'sensory' 'time'},5,subplot(2,3,2:3),traceColor(x,:),angleIdx); %plot regressor angles in current area
    hold(ax,'on');
    
end

% mean standard deviation
subplot(2,3,4);
cMap = arrayShrink(mean(respVar(:,iMod,size(recV,5) + 1,:),4),mask,'split');
cRange = [-max(abs(cMap(:))) max(abs(cMap(:)))];
mapImg = imshow(cMap,cRange);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,colormap_blueblackred);hold on; colorbar
title(['Mean std - ' modLabel{iMod}]);

[ax, ampIdx] = regressorPlot(squeeze(nanmean(respVar(:,iMod,:,:)))',{recLabels{:} 'full' 'motor' 'sensory' 'time'},5,subplot(2,3,5:6));
title(ax,['PSTH std during stimulus/wait period (' modLabel{iMod} ' trials)']);

for x = 1 : size(areaCoord,2)
    
    areas(x) = rectangle('Position',[areaCoord(1,x)-areaCoord(3,x) areaCoord(2,x)-areaCoord(3,x) areaCoord(3,x)*2 areaCoord(3,x)*2],'Curvature',1,'linewidth',2,'linestyle','--','Edgecolor','w','parent',mapImg.Parent);
    text(areaCoord(1,x)-6,areaCoord(2,x),num2str(x),'color','w','FontSize',15,'parent',mapImg.Parent)
    
    traceIdx{x} = find(createCirclesMask(mapImg.CData,(areas(x).Position(1:2) + (areas(x).Position(3)/2)),areas(x).Position(3)/2)); %get mask from current area and find index
    traceIdx{x} = ismember(maskIdx,traceIdx{x});
    
    ax = regressorPlot(squeeze(nanmean(respVar(traceIdx{x},iMod,:,:)))',{recLabels{:} 'full' 'motor' 'sensory' 'time'},5,subplot(2,3,5:6),traceColor(x,:),ampIdx); %plot regressor angles in current area
    hold(ax,'on');
    
end
regressorPlot(squeeze(nanmean(respVar(:,iMod,:,:)))',{recLabels{:} 'full' 'motor' 'sensory' 'time'},5,subplot(2,3,5:6));
