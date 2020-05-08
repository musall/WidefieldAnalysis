%% basic variables
cPath = 'U:\space_managed_data\BpodImager\Animals\'; %data path on the server
trialSegments = [0 60 90 141 171]; %segments to create maps for whole trial regressors
labels = {'Stimulus' 'Wait'}; %segment labels
motorIdx = 16; %index for zero-lag motor regressor
dimCnt = 200; % use only dimCnt dimensions from Vc. This should match what was used in the model

motorLabels = {'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick' 'leverIn' 'fastPupil' 'slowPupil'... %regressors for motor-based reconstruction
    'BaselineMove' 'HandleMove' 'StimulusMove' 'WaitMove' 'piezo' 'whisk' 'bhvVideo'}; 

dataOverview = delayDecRecordings; %get overview for all recordings

% area coordinates. First two values are x/y for center, third entry is radius.
areaCoord(:,1) = [225 520 25]; %coordinates of area 1 (V1).
areaCoord(:,2) = [290 440 20]; %coordinates of area 2 (retrosplinal).
areaCoord(:,3) = [240 350 25]; %coordinates of area 3 (hindlimb).
areaCoord(:,4) = [260 230 25]; %coordinates of area 4 (M2).

%% get reconstructed Vs, modality index, datapath and allOpts
[recV, recLabels, dataPath, allOpts, modIdx] = BpodImager_motorReconstruct('All', 'All', false); %get reconstructed V, used full model
load([dataPath{1} 'interpVc.mat'],'frames')
recV = cat(5,recV{:});

allSnap = NaN(650, 650, length(dataPath), 'single'); %pre-allocate larger data array for brain picture
allU = NaN(650, 650, dimCnt, length(dataPath), 'single'); %pre-allocate larger data array for all Us
allV = NaN(dimCnt, frames, 2, length(dataPath), 'single'); %pre-allocate data array for PSTHs

%% load data
for iAnimals = 1:length(dataPath)
    
    load([dataPath{iAnimals} 'Vc.mat'],'U')
    load([dataPath{iAnimals} 'snapshot_1.mat'])
    U = U(:,:,1:dimCnt);
    
    % re-align U and cut to ensure it fits into larger matrix
    U = Widefield_mapAlign(U,allOpts(iAnimals).opts);
    snap = Widefield_mapAlign(single(snap),allOpts(iAnimals).opts);
    mask = ~isnan(U(:,:,1));
    
    xCut = find(sum(mask,1) > 0);
    xCut = xCut([1 end]);
    xCut = min([xCut(1) size(mask,2)-xCut(end)]);
    
    yCut = find(sum(mask,2) > 0);
    yCut = yCut([1 end]); 
    yCut = min([yCut(1) size(mask,1)-yCut(end)]);
    
    U = U(yCut:end-yCut+1,xCut:end-xCut+1,:);
    snap = snap(yCut:end-yCut+1,xCut:end-xCut+1);
    
    % figure out how to fit into larger array
    yShift = (size(allU,1)-size(U,1))/2;
    xShift = (size(allU,2)-size(U,2))/2;
    
    allSnap(yShift+1:end-yShift,xShift+1:end-xShift,iAnimals) = snap;
    allU(yShift+1:end-yShift,xShift+1:end-xShift,:,iAnimals) = U;
    
    %% get Vc that was used for model and get PSTH for correct vis and aud trials
    load([dataPath{iAnimals} 'interpVc.mat'],'Vc','frames')
    allV(:,:,1,iAnimals) = mean(reshape(Vc(:,modIdx{iAnimals,1}),dimCnt,frames,[]),3);
    allV(:,:,2,iAnimals) = mean(reshape(Vc(:,modIdx{iAnimals,2}),dimCnt,frames,[]),3);
    
end

mask = isnan(nanmean(nanmean(allU,3),4));
allSnap = arrayShrink(allSnap,mask);
allU = arrayShrink(allU,mask);

%% get reconstructed data
fullRec{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
fullRec{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

motorRec{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
motorRec{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

sensoryRec{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
sensoryRec{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

allData{1} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs
allData{2} = NaN(sum(~mask(:)), frames, length(dataPath), 'single'); %pre-allocate data array for PSTHs

% convolve real and reconstructed Vs to create PSTHs
for iAnimals = 1:length(dataPath)

    allData{1}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(allV(:,:,1,iAnimals),4); % original data
    allData{2}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(allV(:,:,2,iAnimals),4); % original data
    
    fullRec{1}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,1,iAnimals,:),5); % full reconstruction
    fullRec{2}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,2,iAnimals,:),5); % full reconstruction
    
    motorRec{1}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,1,iAnimals,ismember(recLabels,motorLabels)),5); % motor reconstruction
    motorRec{2}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,2,iAnimals,ismember(recLabels,motorLabels)),5); % motor reconstruction
    
    sensoryRec{1}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,1,iAnimals,~ismember(recLabels,motorLabels)),5); % non-motor reconstruction
    sensoryRec{2}(:, :, iAnimals) = allU(:,:,iAnimals) * sum(recV(:,:,2,iAnimals,~ismember(recLabels,motorLabels)),5); % non-motor reconstruction

end

%% get predicted R^2 map and cut to size
temp = squeeze(BpodImager_compareModalities('all','fullCorr',[],true));
allCorr = NaN(size(mask,1),size(mask,2),size(temp,3));

yShift = (size(mask,1)-size(temp,1))/2;
xShift = (size(mask,2)-size(temp,2))/2;
if yShift < 0
    yShift = abs(yShift);
    temp = temp(xShift+1:end-xShift,:,:);
    yShift = 0;
end
if xShift < 0
    xShift = abs(xShift);
    temp = temp(:,xShift+1:end-xShift,:);
    xShift = 0;
end
allCorr(yShift+1:end-yShift,xShift+1:end-xShift,:) = temp;
allCorr = arrayShrink(allCorr,mask);

%% show real data, vision and auditory
mapData = nanmean(nanmean(allData{1}(:,trialSegments(3)+1:trialSegments(4),:),3),2);
Widefield_tracePlot(mapData,areaCoord,allData,mask,trialSegments,[-0.005 0.005],colormap_blueblackred); 

% compare real versus fully modeled data
Widefield_tracePlot(nanmean(allCorr,2),areaCoord,{allData{1} fullRec{1}},mask,trialSegments,[0 1],colormap('hot'));

% compare real data, versus motor and 'sensory' component
mapData = nanmean(nanmean(allData{1}(:,trialSegments(3)+1:trialSegments(4),:),3),2);
Widefield_tracePlot(mapData,areaCoord,{allData{1} motorRec{1} sensoryRec{1}},mask,trialSegments,[-0.005 0.005],colormap_blueblackred);
mapData = nanmean(nanmean(allData{2}(:,trialSegments(3)+1:trialSegments(4),:),3),2);
Widefield_tracePlot(mapData,areaCoord,{allData{2} motorRec{2} sensoryRec{2}},mask,trialSegments,[-0.005 0.005],colormap_blueblackred);

% same thing with all data
mapData = nanmean(nanmean(allData{2}(:,trialSegments(3)+1:trialSegments(4),:),3),2);
Widefield_tracePlot(mapData,areaCoord,{cat(3,allData{:}) cat(3,motorRec{:}) cat(3,sensoryRec{:})},mask,trialSegments,[-0.005 0.005],colormap_blueblackred);


%% separate data based on animal expertise
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts

% compare expert vs naive vision
mapData = nanmean(nanmean(allData{1}(:,trialSegments(3)+1:trialSegments(4),visExp),3),2);
Widefield_tracePlot(mapData,areaCoord,{allData{1}(:,:,visExp) allData{1}(:,:,audExp)} ,mask,trialSegments,[-0.005 0.005],colormap_blueblackred); 

% compare expert vs naive audio
mapData = nanmean(nanmean(allData{2}(:,trialSegments(3)+1:trialSegments(4),audExp),3),2);
Widefield_tracePlot(mapData,areaCoord,{allData{2}(:,:,audExp) allData{2}(:,:,visExp)} ,mask,trialSegments,[-0.005 0.005],colormap_blueblackred); 

% compare expert vs naive
mapData = nanmean(nanmean(allData{2}(:,trialSegments(3)+1:trialSegments(4),audExp),3),2);
Widefield_tracePlot(mapData,areaCoord,{cat(3,allData{1}(:,:,visExp),allData{2}(:,:,audExp)) cat(3,allData{2}(:,:,visExp),allData{1}(:,:,audExp))} ,mask,trialSegments,[-0.005 0.005],colormap_blueblackred); 

% compare expert vs naive (motor)
Widefield_tracePlot(mapData,areaCoord,{cat(3,motorRec{1}(:,:,visExp),motorRec{2}(:,:,audExp)) cat(3,motorRec{2}(:,:,visExp),motorRec{1}(:,:,audExp))} ,mask,trialSegments,[-0.005 0.005],colormap_blueblackred); 
Widefield_tracePlot(mapData,areaCoord,{cat(3,sensoryRec{1}(:,:,visExp),sensoryRec{2}(:,:,audExp)) cat(3,sensoryRec{2}(:,:,visExp),sensoryRec{1}(:,:,audExp))} ,mask,trialSegments,[-0.005 0.005],colormap_blueblackred); 



