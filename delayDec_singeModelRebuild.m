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
taskV = BpodImager_motorReconstruct(cPath, 'All', 'task', true, [], true); taskV = taskV{1}; %get reconstructed V, used full model
opMotorV = BpodImager_motorReconstruct(cPath, 'All', 'opMotor', true, [], true); opMotorV = opMotorV{1}; %get reconstructed V, used full model
[spMotorV, recLabels, dataPath, allOpts, modIdx, sideIdx, alignIdx, baseLength, frames, stimTimes] = BpodImager_motorReconstruct(cPath, 'All', 'spMotor', true, [], true);  spMotorV = spMotorV{1}; %get reconstructed V, used full model
dimCnt = size(taskV,1);

%% get allen maps
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

% pre-allocate arrays
allSnap = NaN(sum(~mask(:)), length(dataPath), 'single'); %pre-allocate larger data array for brain picture
allU = NaN(sum(~mask(:)), dimCnt, length(dataPath), 'single'); %pre-allocate larger data array for all Us
allV = NaN(dimCnt, frames, 3, length(dataPath), 'single'); %pre-allocate data array for PSTHs

%% load raw data and get aligned U + averaged Vs for vision, audio and all trials
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
    Vc = bsxfun(@minus, Vc, mean(Vc,2));
    Vc = Vc(:,alignIdx{iAnimals});
    
    allV(:,:,1,iAnimals) = mean(reshape(Vc(:,modIdx{iAnimals,1}),dimCnt,frames,[]),3); %correct vision
    allV(:,:,2,iAnimals) = mean(reshape(Vc(:,modIdx{iAnimals,2}),dimCnt,frames,[]),3); %correct audio
    allV(:,:,3,iAnimals) = mean(reshape(Vc(:,modIdx{iAnimals,3}),dimCnt,frames,[]),3); %all correct
    
end
clear Vc U dimBeta fullR

%% get reconstructed data
cMod = 1; %use all trials
allData = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
opMotorRec = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
spMotorRec = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
taskRec = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

% convolve real and reconstructed Vs to create PSTHs
for iAnimals = 1:length(dataPath)

    allData(:, :, iAnimals) = (allU(:,:,iAnimals) * allV(:,:,cMod,iAnimals)); % original data
    opMotorRec(:, :, iAnimals) = (allU(:,:,iAnimals) * opMotorV(:,:,cMod,iAnimals)); % operant motor
    spMotorRec(:, :, iAnimals) = (allU(:,:,iAnimals) * spMotorV(:,:,cMod,iAnimals)); % spont. motor
    taskRec(:, :, iAnimals) = (allU(:,:,iAnimals) * taskV(:,:,cMod,iAnimals)); % task
    
end

%%
% cMovie = arrayShrink(nanmean(allData,3),mask,'split');
cTitle = {'Real vs. model data' 'All mice - visual trials'};
cRange = [-0.02 0.02];

for iAnimal = 7
% for iAnimal = 1:size(visExp,1)
figure;
set(gcf,'position',[-1920 40 1920 963]);
Cnt = 0;
for iMod = 1:3
    if iMod == 1
        cMovie = arrayShrink(nanmean(taskRec(:,:,iAnimal),3),mask,'split');
        cRange = [-nanmean(abs(cMovie(:))*3) nanmean(abs(cMovie(:)))*3];
    elseif iMod == 2
        cMovie = arrayShrink(nanmean(spMotorRec(:,:,iAnimal),3),mask,'split');
    elseif iMod == 3
        cMovie = arrayShrink(nanmean(opMotorRec(:,:,iAnimal),3),mask,'split');        
    end
        
    cMovie = cMovie - mean(cMovie(:,:,1:15),3);

    for iSegs = 1:5
        Cnt = Cnt + 1;
        subplot(3,5,Cnt)
        mapImg = imshow(nanmean(cMovie(:,:,segIdxRealign{iSegs+1}),3),[cRange]);
        colormap(mapImg.Parent,colormap_blueblackred(256)); axis image
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        title(['RecNr: ' num2str(iAnimal) ' - ' segLabels{iSegs+1}])
        
        hold(mapImg.Parent, 'on');
        for x = 1: length(rightHs)
            plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
        end
        drawnow;

    end
end
% Widefield_tracePlot(nanmean(taskRec(:,85,iAnimal),3),areaCoord,{allData(:,:,iAnimal), taskRec(:,:,iAnimal), spMotorRec(:,:,iAnimal), opMotorRec(:,:,iAnimal)},mask,trialSegments,cRange,colormap_blueblackred,[],['iAnimal: ' num2str(iAnimal)],xRange,[cRange(1)*2 cRange(2)*2], false, true);
Widefield_tracePlot(nanmean(taskRec(:,85,iAnimal),3),areaCoord,{allData(:,:,iAnimal), opMotorRec(:,:,iAnimal), taskRec(:,:,iAnimal), spMotorRec(:,:,iAnimal)},mask,trialSegments,cRange,colormap_blueblackred,[],['iAnimal: ' num2str(iAnimal)],xRange,[cRange(1)*3 cRange(2)*4.5], false, true);
drawnow;
end

% Widefield_tracePlot(nanmean(taskRec(:,85,:),3),areaCoord,{nanmean(allData,3), nanmean(taskRec,3), nanmean(spMotorRec,3), nanmean(opMotorRec,3)},mask,trialSegments,cRange,colormap_blueblackred,[],['iAnimal: ' num2str(iAnimal)],xRange,[-0.01 0.01], false, true);

%%
% cTitle = {'Real vs. model data' 'All mice - visual trials'};
% mapData = arrayShrink(nanmean(taskRec(:, 85, :),3),mask,'split');
% cRange = [-nanmean(abs(allData(:))*2) nanmean(abs(allData(:)))*2];
% Widefield_tracePlot(nanmean(taskRec(:,85,:),3),areaCoord,{allData, fullRec},mask,trialSegments,cRange,colormap_blueblackred,[],cTitle,xRange,[cRange(1) cRange(2)*3.5]);
% 
% cTitle = {'Real vs. model data' 'All mice - visual trials'};
% mapData = arrayShrink(nanmean(taskRec(:, 85, :),3),mask,'split');
% cRange = [-nanmean(abs(allData(:))*2) nanmean(abs(allData(:)))*2];
% Widefield_tracePlot(nanmean(taskRec(:,85,:),3),areaCoord,{allData, spMotorRec},mask,trialSegments,cRange,colormap_blueblackred,[],cTitle,xRange,[cRange(1) cRange(2)*3.5]);
% 
% 
% 
% a = (abs(arrayShrink(nanmean(allData,3),mask,'split')) - abs(arrayShrink(nanmean(spMotorRec,3),mask,'split'))) ./ (abs(arrayShrink(nanmean(allData,3),mask,'split')) + abs(arrayShrink(nanmean(spMotorRec,3),mask,'split')));


