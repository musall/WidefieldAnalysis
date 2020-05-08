function BpodImager_checkGFP(reload)

if ~exist('reload','var') || isempty(reload)
    reload = false;
end

%% select data sets
[dataOverview, ~, sensorLabels, cogLabels, ~, segLabels, ~, cPath] = delayDecRecordings;
dataOverview(:,4) = {'Ai93'};
[dataOverview1, ~, sensorLabels, cogLabels, ~, segLabels, ~, cPath] = delayDecRecordings_GFP;
dataOverview = [dataOverview; dataOverview1];
animals = dataOverview(:,1);

modTypes = {'shCurrent' 'shOther' 'shOtherMotor'};
orgVars = {'fastPupil', 'slowPupil', 'whisk', 'nose', 'face', 'lLick', 'rLick', ...
    'pupils', 'licks', 'bhvVideo', 'Move', 'allMove', 'motorNoVideo', 'motor'}; %make sure to load original version during shOther condition to avoid false results from orthogonalization.
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts

GFPidx = ismember(dataOverview(:,4),'GFP'); %index for GFP controls
CamKidx = ismember(dataOverview(:,4),'gcamp'); %index for gcamp controls
Ai93idx = ismember(dataOverview(:,4),'Ai93'); %index for Ai93 animals

load('allenDorsalMapSM.mat')
areaIdx = {'VISp' 'SSp-ul'}; % V1 - FL
leftHS = ismember(dorsalMaps.sidesSplit,'L'); %index for labels on the left HS
for x = 1:length(areaIdx)
    ind = ismember(dorsalMaps.labelsSplit,areaIdx{x}) & leftHS;
    areaCoord{x} = poly2mask(dorsalMaps.edgeOutlineSplit{ind}(:,2), dorsalMaps.edgeOutlineSplit{ind}(:,1),size(dorsalMaps.hsMask,1),size(dorsalMaps.hsMask,2));
end
leftHS = find(leftHS); %index for labels on the left HS

var1 = 'nose'; %first variable to show unique exp. variance
var2 = 'rGrab'; %second variable to show unique exp. variance
var3 = 'piezo'; %third variable to show unique exp. variance

%% general variables
Paradigm = 'SpatialDisc';
% cPath = 'Y:\data\BpodImager\Animals\'; %Widefield data path on nlsas server
sPath = 'Y:\data\predVariance\'; %local data path to save down results

%% load data
check = true;
if ~reload %check for a saved dataset
    if ~exist(sPath,'dir')
        check = false;
    else
        try
            load([sPath 'fullVarMaps_GFP.mat'],'fullVarMaps', 'meanModel', 'allenMask');
            load([sPath 'fullMaps_GFP.mat'],'fullMaps');
            load([sPath 'noseMaps_GFP.mat'],'noseMaps');
            load([sPath 'handleMaps_GFP.mat'],'handleMaps');
            load([sPath 'pawMaps_GFP.mat'],'pawMaps');
            load([sPath 'expHemo_GFP.mat'],'expHemo');
            load([sPath 'betaTraces_GFP.mat'],'visBeta', 'handleBeta');
            
        catch
            check = false;
        end
    end
    if ~check
        reload = true;
        disp('Could not find saved results. Load predVariance data instead. Maybe go get a coffe in the meantime.')
    end
end

%% go through predVariance folders if reloading that data
% check animal count
animalCnt = size(dataOverview,1);
if reload
    animals = dataOverview(:,1);
    Cnt = 0;
    
    for iAnimals = 1:length(animals)
        
        if Ai93idx(iAnimals)
            cPath = 'U:\smusall\BpodImager\Animals\'; %Widefield data path on grid server
        else
            cPath = 'Y:\data\BpodImager\Animals\'; %Widefield data path on grid server
        end
        
        Cnt = Cnt + 1;
        fPath = [cPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
        disp(fPath); tic;
        
        %load design matrix and check for regressors that were used in the model
        load([fPath 'Snapshot_1.mat']); %load snapshot
        load([fPath 'mask.mat']); %load mask
        allOpts(Cnt) = load([fPath 'opts2.mat']); %load allen alignment file
        
        if Cnt == 1
            % create mask from allen coordinates
            allenMask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
            [x1, y1] = size(allenMask);
            [x2, y2] = size(mask);
            allenMask = allenMask(1:min([x1 x2]), 1:min([y1 y2])); %cut allen mask to size
            redMask = ~dorsalMaps.hsMask(1:min([x1 x2]), 1:min([y1 y2])); %cut reduced mask to size
            [x1, y1] = size(allenMask);
            
            % pre-allocate data matrices
            fullVarMaps = NaN(sum(~allenMask(:)), animalCnt, 'single'); %full data variance map
            fullMaps = NaN(sum(~allenMask(:)), animalCnt, 'single'); %full model map
            noseMaps = NaN(sum(~allenMask(:)), animalCnt, 'single'); %maps for right vis stim
            handleMaps = NaN(sum(~allenMask(:)), animalCnt, 'single'); %maps for right handle grab
            pawMaps = NaN(sum(~allenMask(:)), animalCnt, 'single'); %maps for right handle grab
            expHemo = NaN(1, animalCnt, 'single'); %variance explained by hemo correction
            meanModel = NaN(1, animalCnt, 'single'); %mean model variance
            
        end
        
        % load model data
        load([fPath 'predVariance' filesep 'shCurrent' filesep 'fullcorr.mat'], 'cMap'); %load current maps
        cMap = arrayShrink(cMap.^2,mask,'split');
        cMap = alignAllenTransIm(cMap,allOpts(Cnt).opts.transParams);
        fullMaps(:, Cnt) = arrayShrink(cMap(1:x1,1:y1), allenMask, 'merge');
        cMap = arrayShrink(cMap(1:x1,1:y1),redMask);
        meanModel(1,Cnt) = nanmean(cMap(:));
        
        % load nose map
        load([fPath 'predVariance' filesep 'shCurrent' filesep var1 'corr.mat'], 'cMap'); %load current visual map
        cMap = arrayShrink(cMap.^2,mask,'split');
        cMap = alignAllenTransIm(cMap,allOpts(Cnt).opts.transParams);
        noseMaps(:, Cnt) = fullMaps(:, Cnt) - arrayShrink(cMap(1:x1,1:y1), allenMask, 'merge');

        % load handle map
        load([fPath 'predVariance' filesep 'shCurrent' filesep var2 'corr.mat'], 'cMap'); %load current visual map
        cMap = arrayShrink(cMap.^2,mask,'split');
        cMap = alignAllenTransIm(cMap,allOpts(Cnt).opts.transParams);
        handleMaps(:, Cnt) = fullMaps(:, Cnt) - arrayShrink(cMap(1:x1,1:y1), allenMask, 'merge');
             
        % load paw map
        load([fPath 'predVariance' filesep 'shCurrent' filesep var3 'corr.mat'], 'cMap'); %load current visual map
        cMap = arrayShrink(cMap.^2,mask,'split');
        cMap = alignAllenTransIm(cMap,allOpts(Cnt).opts.transParams);
        pawMaps(:, Cnt) = fullMaps(:, Cnt) - arrayShrink(cMap(1:x1,1:y1), allenMask, 'merge');
            
        % load beta maps
        load([fPath 'Vc.mat'], 'U');
        U = alignAllenTransIm(U,allOpts(Cnt).opts.transParams);
        U = U(1:size(allenMask,1), 1:size(allenMask,2), :); 
        
        load([fPath 'dimBeta.mat'], 'dimBeta'); %load current maps
        load([fPath 'regData.mat'], 'idx' ,'trialIdx', 'recIdx', 'recLabels','gaussShift');
        
        temp = dimBeta(recIdx(~idx) == find(ismember(recLabels, var1)), :);
        visBeta(1:size(temp,1), Cnt) = nanmean(arrayShrink(U, ~areaCoord{1}, 'merge')  * temp');
        
        temp = dimBeta(recIdx(~idx) == find(ismember(recLabels, var2)), :);
        handleBeta(1:size(temp,1), Cnt) = nanmean(arrayShrink(U, ~areaCoord{2}, 'merge')  * temp');        
                
        % load explaiend hemo var
        load([fPath 'HemoCorrection.mat'], 'hemoVar'); %load current visual map
        expHemo(Cnt) = hemoVar;
        
        % load explaiend hemo var
        load([fPath 'HemoCorrection.mat'], 'hemoVar'); %load current visual map
        
        % load data and compute its variance
        load([fPath 'Vc.mat'], 'U', 'Vc'); %load current data
        Vc = reshape(Vc,size(Vc,1),[]);
        U = alignAllenTransIm(U,allOpts(Cnt).opts.transParams);
        U = arrayShrink(U(1:x1,1:y1,:), allenMask);
        covVc = cov(Vc');  % S x S
        fullVarMaps(:, Cnt) = sum((U * covVc) .* U, 2)';  % 1 x P
        toc;
    end
    clear cMap Vc U
    
    %% save down results
    if ~exist(sPath,'dir')
        mkdir(sPath);
    end
    save([sPath 'fullVarMaps_GFP.mat'],'fullVarMaps', 'meanModel', 'allenMask');
    save([sPath 'fullMaps_GFP.mat'],'fullMaps');
    save([sPath 'noseMaps_GFP.mat'],'noseMaps');
    save([sPath 'handleMaps_GFP.mat'],'handleMaps');
    save([sPath 'pawMaps_GFP.mat'],'pawMaps');
    save([sPath 'expHemo_GFP.mat'],'expHemo');
    save([sPath 'betaTraces_GFP.mat'],'visBeta', 'handleBeta');
            
end

%% average hemo correction performance
h = figure;
% errorbar([mean(expHemo(GFPidx)) mean(expHemo(Ai93idx)) mean(expHemo(CamKidx))],[sem(expHemo(GFPidx)) sem(expHemo(Ai93idx)) sem(expHemo(CamKidx))], 'k.', 'linewidth',2)
% hold; bar([mean(expHemo(GFPidx)) mean(expHemo(Ai93idx)) mean(expHemo(CamKidx))]);

temp = NaN(size(expHemo,2), 3);
temp(GFPidx, 1) = expHemo(GFPidx);
temp(Ai93idx, 2) = expHemo(Ai93idx);
temp(CamKidx, 3) = expHemo(CamKidx);

boxplot(temp);
axis square; 
set(h.Children,'xTick',1:3)
set(h.Children,'xTickLabel',{'GFP' 'Ai93' 'CamK'})
ylabel('Avg. explained variance');
title('Hemodynamic correction performance'); ylim([0 100]);

%% full model plot
figure;
maxRange = 0.5;
% maxRange = prctile(nanmean(fullMaps(:, Ai93idx),2),99);

subplot(1,3,1);
mapImg = imshow(arrayShrink(squeeze(nanmean(fullMaps(:, GFPidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title(['GFP mice - full model'])
hold(mapImg.Parent, 'on');
for x = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{(x)}(:,2), dorsalMaps.edgeOutlineSplit{(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,3,2);
mapImg = imshow(arrayShrink(squeeze(nanmean(fullMaps(:, Ai93idx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title(['GFP mice - full model'])
hold(mapImg.Parent, 'on');
for x = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{(x)}(:,2), dorsalMaps.edgeOutlineSplit{(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,3,3);
mapImg = imshow(arrayShrink(squeeze(nanmean(fullMaps(:, CamKidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title(['GFP mice - full model'])
hold(mapImg.Parent, 'on');
for x = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{(x)}(:,2), dorsalMaps.edgeOutlineSplit{(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end
    

%% full model average R2
h = figure;
% errorbar([mean(meanModel(GFPidx)) mean(meanModel(Ai93idx)) mean(meanModel(CamKidx))],[sem(meanModel(GFPidx)) sem(meanModel(Ai93idx)) sem(meanModel(CamKidx))], 'k.', 'linewidth',2)
% hold; bar([mean(meanModel(GFPidx)) mean(meanModel(Ai93idx)) mean(meanModel(CamKidx))]);
temp = NaN(size(meanModel,2), 3);
temp(GFPidx, 1) = meanModel(GFPidx);
temp(Ai93idx, 2) = meanModel(Ai93idx);
temp(CamKidx, 3) = meanModel(CamKidx);

boxplot(temp);
axis square; 
set(h.Children,'xTick',1:3)
set(h.Children,'xTickLabel',{'GFP' 'Ai93' 'CamK'})
ylabel('Avg. explained variance');
title('Full model comparison')

%% full variance plot
figure;
% maxRange = 0.5;
maxRange = prctile(nanmean(fullVarMaps(:, Ai93idx),2),99);

subplot(1,3,1);
mapImg = imshow(arrayShrink(squeeze(nanmean(fullVarMaps(:, GFPidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title(['GFP mice - all  variance'])

subplot(1,3,2);
mapImg = imshow(arrayShrink(squeeze(nanmean(fullVarMaps(:, Ai93idx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title(['GFP mice - all  variance'])

subplot(1,3,3);
mapImg = imshow(arrayShrink(squeeze(nanmean(fullVarMaps(:, CamKidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title(['GFP mice - all  variance'])


%% full model explained variance plot
figure;
% maxRange = 0.5;
maxRange = prctile(nanmean(fullMaps(:, Ai93idx).*fullVarMaps(:, Ai93idx),2),99);

subplot(1,3,1);
mapImg = imshow(arrayShrink(squeeze(nanmean(fullMaps(:, GFPidx).*fullVarMaps(:, GFPidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title(['GFP mice - absolute explained  variance'])

subplot(1,3,2);
mapImg = imshow(arrayShrink(squeeze(nanmean(fullMaps(:, Ai93idx).*fullVarMaps(:, Ai93idx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title(['GFP mice - absolute explained  variance'])

subplot(1,3,3);
mapImg = imshow(arrayShrink(squeeze(nanmean(fullMaps(:, CamKidx).*fullVarMaps(:, CamKidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title(['GFP mice - absolute explained  variance'])


%% nose unique variance
figure;
maxRange = 0.03;
% maxRange = prctile(nanmean(noseMaps(:, Ai93idx),2),99);

subplot(1,3,1);
mapImg = imshow(arrayShrink(squeeze(nanmean(noseMaps(:, GFPidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title('GFP mice - nose exp. variance')

subplot(1,3,2);
mapImg = imshow(arrayShrink(squeeze(nanmean(noseMaps(:, Ai93idx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title('Ai93 mice - nose exp. variance')

subplot(1,3,3);
mapImg = imshow(arrayShrink(squeeze(nanmean(noseMaps(:, CamKidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title('CamK mice - nose exp. variance')

%% handle unique variance
figure;
maxRange = 0.02;
% maxRange = prctile(nanmean(handleMaps(:, Ai93idx),2),99);

subplot(1,3,1);
mapImg = imshow(arrayShrink(squeeze(nanmean(handleMaps(:, GFPidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title('GFP mice - handle exp. variance')
hold(mapImg.Parent, 'on');
for x = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{(x)}(:,2), dorsalMaps.edgeOutlineSplit{(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,3,2);
mapImg = imshow(arrayShrink(squeeze(nanmean(handleMaps(:, Ai93idx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title('Ai93 mice - handle exp. variance')
hold(mapImg.Parent, 'on');
for x = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{(x)}(:,2), dorsalMaps.edgeOutlineSplit{(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,3,3);
mapImg = imshow(arrayShrink(squeeze(nanmean(handleMaps(:, CamKidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title('CamK mice - handle exp. variance')
hold(mapImg.Parent, 'on');
for x = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{(x)}(:,2), dorsalMaps.edgeOutlineSplit{(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

%% hindpaw unique variance
figure;
maxRange = 0.01;
% maxRange = prctile(nanmean(handleMaps(:, Ai93idx),2),99);

subplot(1,3,1);
mapImg = imshow(arrayShrink(squeeze(nanmean(pawMaps(:, GFPidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title('GFP mice - hindpaw exp. variance')
hold(mapImg.Parent, 'on');
for x = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{(x)}(:,2), dorsalMaps.edgeOutlineSplit{(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,3,2);
mapImg = imshow(arrayShrink(squeeze(nanmean(pawMaps(:, Ai93idx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title('Ai93 mice - hindpaw exp. variance')
hold(mapImg.Parent, 'on');
for x = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{(x)}(:,2), dorsalMaps.edgeOutlineSplit{(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

subplot(1,3,3);
mapImg = imshow(arrayShrink(squeeze(nanmean(pawMaps(:, CamKidx),2)), allenMask, 'split'),[0 maxRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); axis image
title('CamK mice - hindpaw exp. variance')
hold(mapImg.Parent, 'on');
for x = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{(x)}(:,2), dorsalMaps.edgeOutlineSplit{(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
end

%% beta traces
figure;
subplot(1,2,1); hold on;
stdshade((visBeta(2:45, GFPidx) - visBeta(1, GFPidx))', 'b', (1:44)*(1/30),  0.5);
stdshade((visBeta(2:45, Ai93idx) - visBeta(1, Ai93idx))', 'k', (1:44)*(1/30), 0.5);
stdshade((visBeta(2:45, CamKidx) - visBeta(1, CamKidx))', 'r', (1:44)*(1/30), 0.5);
legend({'GFP' 'Ai93' 'CamK'}); axis square;
xlabel('Time (s)'); ylabel('dF/F');

subplot(1,2,2); hold on;
stdshade((handleBeta(2:45, GFPidx) - handleBeta(1, GFPidx))', 'b', (1:44)*(1/30), 0.5);
stdshade((handleBeta(2:45, Ai93idx) - handleBeta(1, Ai93idx))', 'k', (1:44)*(1/30), 0.5);
stdshade((handleBeta(2:45, CamKidx) - handleBeta(1, CamKidx))', 'r', (1:44)*(1/30), 0.5);
legend({'GFP' 'Ai93' 'CamK'}); axis square
xlabel('Time (s)'); ylabel('dF/F');

