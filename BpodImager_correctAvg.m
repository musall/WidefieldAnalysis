%% basic variables
cPath = 'U:\space_managed_data\BpodImager\Animals\'; %data path on the server
trialSegments = [0 55 80 130 160]; %segments to create maps for whole trial regressors
motorIdx = 16; %index for zero-lag motor regressor
dimCnt = 200; % use only dimCnt dimensions from Vc. This should match what was used in the model
baseRange = 1:15; %baseline frames
areaIdx = {'VISp' 'RSPd' 'SSp-ll' 'MOs'}; % V1 - RS - HL - M2
normData = false; %flag to zscore data
sRate = 30; %sampling rate in Hz
runModel = true; %run model in motor-subtracted data

%% session data
[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, segIdxRealign] = delayDecRecordings;
sensorLabels = [sensorLabels cogLabels]; %combine sensory and cognitive labels, so sensory is really 'non-motor' or 'task' here
segIdxRealign{2} = 46:65;
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts

%% get indices for modalities, datapath and allOpts
[dataPath, allOpts, modIdx, sideIdx, alignIdx, baseLength, frames, stimTimes] = BpodImager_avgIndex('All', false); %get average Vc after removing motor model

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
xRange = [(-baseLength : 1 : 0)./sRate, (1 : frames-baseLength-1) ./sRate]; %time vector for x-axis. Stimulus onset is at 0s.
allSnap = NaN(sum(~mask(:)), length(dataPath), 'single'); %pre-allocate larger data array for brain picture
allU = NaN(sum(~mask(:)), dimCnt, length(dataPath), 'single'); %pre-allocate larger data array for all Us

allData{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs.
allData{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
allData{3} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

recData{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for task PSTHs.
recData{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for task PSTHs
recData{3} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for task PSTHs

dpData{1} = NaN(sum(~mask(:)), 2, length(dataPath)); %dprime maps for expertise. original data.
dpData{2} = NaN(sum(~mask(:)), 2, length(dataPath)); %dprime maps for expertise. motor-subtract data.

%% load U and snapshots for each recording and reconstruct task-related data
recV = NaN(dimCnt, frames, 3, length(dataPath), 2);
for iAnimals = 1:length(dataPath)
    
    % get Vc and Vm. Align to create PSTHs and correct for baseline.
    load([dataPath{iAnimals} 'interpVc.mat'],'Vc')
    load([dataPath{iAnimals} 'interpVmotor.mat'],'Vm')

    Vc = Vc(:, alignIdx{iAnimals}); %reduce Vc to only include aligned baseline and poststim data
    Vm = Vm(:, alignIdx{iAnimals}); %reduce Vm to only include aligned baseline and poststim data
    Vm = Vc - Vm; %remove Vmotor from Vc
    
    %% get Vc baseline and correct both datasets accordingly
    Vc = reshape(Vc,size(Vc,1),frames,[]);
    temp = squeeze(mean(mean(Vc(:,1:sRate/2,:),3),2));
    Vc = reshape(Vc, size(Vc,1), []);
    Vc = bsxfun(@minus, Vc, temp);
    
    Vm = reshape(Vm,size(Vm,1),frames,[]);
    temp = squeeze(mean(mean(Vm(:,1:sRate/2,:),3),2));
    Vm = reshape(Vm, size(Vm,1), []);
    Vm = bsxfun(@minus, Vm, temp);

    for iMod = 1:3 %different modalities, 1 is corr. vision, 2 is corr. audio, 3 is all corr. trials. 4:6 is the same thing but for all trials.
        temp = reshape(Vc(:, modIdx{iAnimals,iMod}),size(Vc,1),frames,[]); %get psth
        temp = mean(temp,3);
        recV(:,:,iMod,iAnimals,1) = temp; %PSTH from original data
        
        temp = reshape(Vm(:,modIdx{iAnimals,iMod}),size(Vm,1),frames,[]); %get psth
        temp = mean(temp,3);
        recV(:,:,iMod,iAnimals,2) = temp; %PSTH from motor-removed data
    end
    
    % get U and do spatial reconstruction.
    load([dataPath{iAnimals} 'Vc.mat'],'U')
    load([dataPath{iAnimals} 'snapshot_1.mat'])
    U = U(:,:,1:size(recV,1));
    
    snap = alignAllenTransIm(single(snap),allOpts(iAnimals).opts.transParams);
    allSnap(:,iAnimals) = arrayShrink(snap(1:size(mask,1),1:size(mask,2)),mask);
    U = alignAllenTransIm(U,allOpts(iAnimals).opts.transParams);
    U = arrayShrink(U(1:size(mask,1),1:size(mask,2),:),mask);
          
    if visExp(iAnimals)
        cInd = [1 2];
    else
        cInd = [2 1];
    end
    
    allData{1}(:, :, iAnimals) = U * recV(:,:,3,iAnimals,1); % all correct PSTH
    allData{2}(:, :, iAnimals) = U * recV(:,:,cInd(1),iAnimals,1); % expert PSTH
    allData{3}(:, :, iAnimals) = U * recV(:,:,cInd(2),iAnimals,1); % novice PSTH
    
    recData{1}(:, :, iAnimals) = U * recV(:,:,3,iAnimals,2); % all correct PSTH
    recData{2}(:, :, iAnimals) = U * recV(:,:,cInd(1),iAnimals,2); % expert PSTH
    recData{3}(:, :, iAnimals) = U * recV(:,:,cInd(2),iAnimals,2); % novice PSTH
    
    %% compute standard deviation and dPrime maps
    x = [3 5]; % segments from segIdx for d' calculation
    cIdx = unique(ceil(find(modIdx{iAnimals,3})/frames)); %index for all unisensory, succesful trial
    
    Vc = reshape(Vc,size(Vc,1),frames,[]);
    Vm = reshape(Vm,size(Vm,1),frames,[]);

    for iSegs = 1:2
        % get dprime for original data
        cData = squeeze(mean(Vc(:, segIdxRealign{x(iSegs)},cIdx),2));
        covV = cov(cData');  % S x S
        varP = sum((U * covV) .* U, 2);  % 1 x P
        stdP = sqrt(varP); %standard deviation map for current segment in all correct trials
        
        expResponse = mean(allData{2}(:, segIdxRealign{x(iSegs)}, iAnimals),2); %response to expert modality
        naiveResponse = mean(allData{3}(:, segIdxRealign{x(iSegs)}, iAnimals),2); %response to naive modality
        dpData{1}(:,iSegs,iAnimals) = (expResponse - naiveResponse) ./ stdP; % d' original data. expertise.
        
        % get dprime for motor-subtrated data
        cData = squeeze(mean(Vm(:, segIdxRealign{x(iSegs)},cIdx),2));
        covV = cov(cData');  % S x S
        varP = sum((U * covV) .* U, 2);  % 1 x P
        stdP = sqrt(varP); %standard deviation map for current segment in all correct trials
        
        expResponse = mean(recData{2}(:, segIdxRealign{x(iSegs)}, iAnimals),2); %response to expert modality
        naiveResponse = mean(recData{3}(:, segIdxRealign{x(iSegs)}, iAnimals),2); %response to naive modality
        dpData{2}(:,iSegs,iAnimals) = (expResponse - naiveResponse) ./ stdP; % d' original data. expertise.
    end
    if rem(iAnimals, round(length(dataPath)/5)) == 0
        fprintf(1, 'Done. Current recording is %d out of %d\n', iAnimals, length(dataPath));
    end
end

%%
figure;
cRange = [0 0.0075];
for iSegs = 1:5
    subplot(2,5,iSegs); hold off;
    cMovie = (arrayShrink(nanmean(abs(allData{1}(:,segIdxRealign{iSegs+1},:)),3),mask,'split'));
    cMap = nanmean((cMovie),3);
    mapImg = imshow(cMap,cRange);
    colormap(mapImg.Parent,'jet'); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    hold(mapImg.Parent, 'on');
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end

cRange = [0 0.0015];
for iSegs = 1:5
    subplot(2,5,5+iSegs); hold off;
    cMovie = (arrayShrink(nanmean(abs(recData{1}(:,segIdxRealign{iSegs+1},:)),3),mask,'split'));
    cMap = nanmean((cMovie),3);
    mapImg = imshow(cMap,cRange);
    colormap(mapImg.Parent,'jet'); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    hold(mapImg.Parent, 'on');
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end


%% show motor-removed data
% cMovie = nanmean(arrayShrink(cat(3, recData{2}(:,:,visExp), recData{3}(:,:,audExp)),mask,'split'),4);
cMovie = (arrayShrink(nanmean(recData{1},3),mask,'split'));
cRange = [-nanmean(abs(cMovie(:))*3) nanmean(abs(cMovie(:)))*3];
figure;
for iSegs = 1:5
    subplot(2,5,iSegs)
    mapImg = imshow(nanmean(cMovie(:,:,segIdxRealign{iSegs+1}),3),[cRange]);
    colormap(mapImg.Parent,colormap_blueblackred); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['All mice. Non-motor rebuild: ' segLabels{iSegs+1}])
    
    hold(mapImg.Parent, 'on');
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
end


%%
figure
Cnt = 0;
for xx = 1:2
    for y = 1:2
        Cnt = Cnt + 1;
        subplot(2,2,Cnt)
        mapImg = imshow((arrayShrink(nanmean(dpData{xx}(:,y,:),3),mask,'split')),[-0.25 0.25]);
        colormap(mapImg.Parent,'jet'); axis image
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        
%         hold(mapImg.Parent, 'on');
%         for x = 1: length(rightHs)
%             plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
%         end
    end
end



