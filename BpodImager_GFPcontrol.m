function BpodImager_GFPcontrol(reload)

if ~exist('reload','var') || isempty(reload)
    reload = false;
end

%% select data sets
controlType = 'gcamp';
[dataOverview, motorLabels, sensorLabels, cogLabels, ~, segLabels, ~, cPath] = delayDecRecordings_GFP;
dataOverview = dataOverview(ismember(dataOverview(:,4), controlType), :);
animals = dataOverview(:,1);

visExp = ismember(dataOverview(1:size(dataOverview,1),2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(1:size(dataOverview,1),2),'Audio'); %index for auditory experts
% modTypes = {'shCurrent' 'shOther' 'shTaskCurrent' 'shOtherMotor'};
modTypes = {'shCurrent' 'shOther' 'shOtherMotor'};
orgVars = {'fastPupil', 'slowPupil', 'whisk', 'nose', 'face', 'lLick', 'rLick', ...
    'pupils', 'licks', 'bhvVideo', 'Move', 'allMove', 'motorNoVideo', 'motor'}; %make sure to load original version during shOther condition to avoid false results from orthogonalization.
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts

%% general variables
Paradigm = 'SpatialDisc';
% cPath = 'U:\smusall\BpodImager\Animals\'; %Widefield data path on grid server
% cPath = 'Y:\data\BpodImager\Animals\'; %Widefield data path on nlsas server
sPath = 'M:\BpodImager\predVariance\'; %local data path to save down results

%% load data
check = true;
if ~reload %check for a saved dataset
    if ~exist(sPath,'dir')
        check = false;
    else
        try
            load([sPath controlType '_corrMaps_GFP.mat'],'corrMaps','allenMask','dorsalMaps');
            load([sPath controlType '_meanRsq_GFP.mat'],'aRsq','aTimeRsq','nRecLabels');
            load([sPath controlType '_segMovies_GFP.mat'],'segMovies');
            load([sPath controlType '_fullMovies_GFP.mat'],'fullMovies','taskMovies','motorMovies');
            load([sPath controlType '_fullCorrMaps_GFP.mat'],'fullCorrMaps');
            load([sPath controlType '_fullSegMovies_GFP.mat'],'fullSegMovies');
            load([sPath controlType '_otherMsegMovies_GFP.mat'],'otherMsegMovies');
            load([sPath controlType '_otherMcorrMaps_GFP.mat'],'otherMcorrMaps', 'oMotorLabels', 'otherMotorLabels');
            
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
        regCnt = 0;
        
        Cnt = Cnt + 1;
        fPath = [cPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
        disp(fPath); tic;
        
        %load design matrix and check for regressors that were used in the model
        load([fPath 'regData.mat'],'recIdx','idx','recLabels')
        usedR{Cnt} = recLabels(unique(recIdx(~idx))); %regressors that were used in the model
        
        load([fPath 'Snapshot_1.mat']); %load snapshot
        load([fPath 'mask.mat']); %load mask
        allOpts(Cnt) = load([fPath 'opts2.mat']); %load allen alignment file
        
        if Cnt == 1
            % create mask from allen coordinates
            load('allenDorsalMapSM.mat')
            allenMask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
            [x1, y1] = size(allenMask);
            [x2, y2] = size(mask);
            allenMask = allenMask(1:min([x1 x2]), 1:min([y1 y2])); %cut allen mask to size
            redMask = ~dorsalMaps.hsMask(1:min([x1 x2]), 1:min([y1 y2])); %cut reduced mask to size
%             redMask = allenMask; %cut reduced mask to size
            [x1, y1] = size(allenMask);
            
            % get labels for all regressors
            load([fPath 'predVariance' filesep modTypes{1} filesep 'fullcorr.mat'],'recLabels', 'segMovie','cMovie'); %load current labels
            load([fPath 'predVariance' filesep 'extraGroups.mat'],'extraGroups'); %load extra group labels
            load([fPath 'predVariance' filesep 'oMotorLabels.mat'],'oMotorLabels'); %load other motor labels
            nRecLabels = [recLabels{:} extraGroups(1,:)]; %labels of all regressors
            nRecLabels = nRecLabels(~strcmpi(nRecLabels,'leverIn')); %don't use leverIn regressor
            
            % make recLabels where modality is swapped for expertise
            altRecLabels = nRecLabels;
            altRecLabels(ismember(nRecLabels,'visReward')) = {'expReward'};
            altRecLabels(ismember(nRecLabels,'audReward')) = {'novReward'};
            
            % pre-allocate data matrices
            fullCorrMaps = NaN(sum(~allenMask(:)), animalCnt, 'single');
            fullSegMovies = NaN(sum(~allenMask(:)), animalCnt, size(segMovie,2), 'single');
            fullMovies = NaN(sum(~allenMask(:)), animalCnt, size(cMovie,2), 'single');
            taskMovies = NaN(sum(~allenMask(:)), animalCnt, size(cMovie,2), 'single');
            motorMovies = NaN(sum(~allenMask(:)), animalCnt, size(cMovie,2), 'single');
            
            corrMaps = NaN(sum(~allenMask(:)), animalCnt, length(nRecLabels), length(modTypes)-1, 'single');
            segMovies = NaN(sum(~allenMask(:)), animalCnt, length(nRecLabels), size(segMovie,2), length(modTypes)-1, 'single');
            
            otherMcorrMaps = NaN(sum(~allenMask(:)), animalCnt, length(oMotorLabels), 'single');
            otherMsegMovies = NaN(sum(~allenMask(:)), animalCnt, length(oMotorLabels), size(segMovie,2), 'single');
            
            aRsq = NaN(length(nRecLabels) + 1, animalCnt, length(modTypes), 'single');
            aTimeRsq = NaN(length(nRecLabels) + 1, size(segMovie,2), animalCnt, length(modTypes), 'single');
        end
        
        
        for iRuns = 1:length(modTypes)
            
            if iRuns == 1 %load full model data only once - running through different versions makes no difference here
                
                load([fPath 'predVariance' filesep modTypes{iRuns} filesep 'fullcorr.mat'],'cMap','segMovie','cMovie','recLabels'); %load current data
                
                % align cMap to allen coordinates
                cMap = arrayShrink(cMap.^2,mask,'split');
                cMap = alignAllenTransIm(cMap,allOpts(Cnt).opts.transParams);
                
                redMap = arrayShrink(cMap(1:x1,1:y1),redMask);
                aRsq(1,Cnt,1) = nanmean(redMap(:));

                cMap = arrayShrink(cMap(1:x1,1:y1),allenMask);
                fullCorrMaps(:, Cnt) = cMap;
                
                % align segMovie to allen coordinates
                segMovie = arrayShrink(segMovie.^2,mask,'split');
                segMovie = alignAllenTransIm(segMovie,allOpts(Cnt).opts.transParams);
                
                redMovie = arrayShrink(segMovie(1:x1,1:y1,:),redMask);
                aTimeRsq(1,:,Cnt,1) = nanmean(redMovie);
                
                segMovie = arrayShrink(segMovie(1:x1,1:y1,:),allenMask);
                fullSegMovies(:, Cnt, :) = segMovie;
                
                % align cMovie to allen coordinates
                cMovie = arrayShrink(cMovie.^2,mask,'split');
                cMovie = alignAllenTransIm(cMovie,allOpts(Cnt).opts.transParams);
                cMovie = arrayShrink(cMovie(1:x1,1:y1,:),allenMask);
                fullMovies(:, Cnt, :) = cMovie;
                
            end
            
            for iRegs = 1:length(nRecLabels)
                try
                    %check if current regressor is task or motor model and load cMovie for that case
                    if (strcmpi(modTypes{iRuns},'shCurrent') && strcmpi(nRecLabels{iRegs},'motor')) || ... %task model. load movie
                            (strcmpi(modTypes{iRuns},'shOther') && strcmpi(nRecLabels{iRegs},'motor')) %motor model. load movie.
                        
                        load([fPath 'predVariance' filesep modTypes{iRuns} filesep nRecLabels{iRegs} 'corr.mat'],'cMap','segMovie','cMovie'); %load current data
                        cMovie = arrayShrink(cMovie.^2,mask,'split');
                        cMovie = alignAllenTransIm(cMovie,allOpts(Cnt).opts.transParams);
                        cMovie = arrayShrink(cMovie(1:x1,1:y1,:),allenMask);
                        if strcmpi(modTypes{iRuns},'shCurrent')
                            taskMovies(:, Cnt, :) = cMovie;
                        else
                            motorMovies(:, Cnt, :) = cMovie;
                        end
                    else
                        if ~ismember(nRecLabels{iRegs}, [usedR{Cnt} extraGroups(1,:)]) && strcmpi(modTypes{iRuns},'shOther') % if regressor was not used and other regressors were shuffled, set all results to zero
                            load([fPath 'predVariance' filesep 'shCurrent' filesep 'fullcorr.mat'],'cMap','segMovie'); %load full model instead of current regressor data
                            cMap = cMap - cMap; %set to 0
                            segMovie = segMovie - segMovie; %set to 0
                        elseif (ismember(nRecLabels{iRegs}, orgVars) && strcmpi(modTypes{iRuns},'shOther')) || strcmpi(modTypes{iRuns},'shOtherMotor') % use original for standalone video regressors
                            load([fPath 'predVariance' filesep modTypes{iRuns} filesep 'org' nRecLabels{iRegs} 'corr.mat'],'cMap','segMovie'); %load current data
                        else
                            load([fPath 'predVariance' filesep modTypes{iRuns} filesep nRecLabels{iRegs} 'corr.mat'],'cMap','segMovie'); %load current data
                        end
                    end
                    
                    segMovie = arrayShrink(segMovie.^2,mask,'split');
                    segMovie = alignAllenTransIm(segMovie,allOpts(Cnt).opts.transParams);
                    segMovie = arrayShrink(segMovie(1:x1,1:y1,:),allenMask);
                    
                    cMap = arrayShrink(cMap.^2,mask,'split');
                    cMap = alignAllenTransIm(cMap,allOpts(Cnt).opts.transParams);
                    redMap = arrayShrink(cMap(1:x1,1:y1),redMask);
                    cMap = arrayShrink(cMap(1:x1,1:y1),allenMask);
                    
                    % change current regressor index if modality should point to animal experitise instead.
                    cInd = iRegs;
                    if strcmpi(nRecLabels{iRegs},'visReward') || strcmpi(nRecLabels{iRegs},'audReward')
                        if (visExp(iAnimals) && strcmpi(nRecLabels{iRegs},'visReward')) || (audExp(iAnimals) && strcmpi(nRecLabels{iRegs},'audReward'))
                            cInd = find(ismember(altRecLabels, 'expReward')); %change index to point to expert instead of modality
                        else
                            cInd = find(ismember(altRecLabels, 'novReward')); %change index to point to novice instead of modality
                        end
                    end
                    aTimeRsq(cInd + 1, :, Cnt, iRuns) = nanmean(segMovie);
                    aRsq(cInd + 1, Cnt, iRuns) = nanmean(redMap(:));
                    
                    if ~strcmpi(modTypes{iRuns},'shOtherMotor')
                        segMovies(:, Cnt, cInd, :, iRuns) = segMovie;
                        corrMaps(:, Cnt, cInd, iRuns) = cMap;
                    elseif any(ismember(oMotorLabels,nRecLabels{iRegs}))
                        regCnt = regCnt + 1;
                        otherMsegMovies(:, Cnt, regCnt, :) = segMovie;
                        otherMcorrMaps(:, Cnt, regCnt) = cMap;
                        otherMotorLabels{regCnt} = nRecLabels{iRegs}; %this is a control. should equal oMotorLabels in the end.
                    end
                end
            end
        end
        toc;
    end
    clear segMovie cMap
    
    %% save down results
    if ~exist(sPath,'dir')
        mkdir(sPath);
    end
    save([sPath controlType '_fullCorrMaps_GFP.mat'],'fullCorrMaps');
    save([sPath controlType '_fullSegMovies_GFP.mat'],'fullSegMovies');
    save([sPath controlType '_fullMovies_GFP.mat'],'fullMovies','taskMovies','motorMovies', '-v7.3');
    save([sPath controlType '_corrMaps_GFP.mat'],'corrMaps','allenMask','dorsalMaps');
    save([sPath controlType '_segMovies_GFP.mat'],'segMovies', '-v7.3');
    save([sPath controlType '_otherMsegMovies_GFP.mat'],'otherMsegMovies', '-v7.3');
    save([sPath controlType '_otherMcorrMaps_GFP.mat'],'otherMcorrMaps', 'oMotorLabels', 'otherMotorLabels');
    save([sPath controlType '_meanRsq_GFP.mat'],'aRsq','aTimeRsq','nRecLabels');
end

%% check if shOtherMotor was loaded correctly
if length(oMotorLabels) ~= sum(ismember(oMotorLabels,otherMotorLabels))
    error('oMotorLabels and otherMotorLabels are not the same. Problem with loading shOtherMotor variables?')
end

%%
segLabels = [{'all'} segLabels];
% fullRecLabels = [{'full'} nRecLabels]; %add full moodel label
altRecLabels = nRecLabels;
altRecLabels(ismember(nRecLabels,'visReward')) = {'expReward'};
altRecLabels(ismember(nRecLabels,'audReward')) = {'novReward'};
fullRecLabels = [{'full'} altRecLabels]; %add full moodel label

%% crossvalidated R2 maps - full model
cMap = arrayShrink(nanmean(fullCorrMaps,2),allenMask,'split');
cMovie = arrayShrink(nanmean(fullSegMovies,2),allenMask,'split');
cMovie = cat(3,cMap,squeeze(cMovie)); %combine full correlation map with different segments
% cRange = max(abs(cMap(:))).*0.75;
cRange = 0.4;

figure;
for iSegs = 1:size(cMovie,3)
    subplot(2,4,iSegs)
    mapImg = imshow(cMovie(:,:,iSegs),[0 cRange]);axis image;
    set(mapImg,'AlphaData',~allenMask); %make NaNs transparent.
    colormap(mapImg.Parent,inferno(256));
    title(['Full model, R^2 - ' segLabels{iSegs}])
end

%% make overview figure for overall task- and motor explained variance
fullMap = arrayShrink(nanmean(fullCorrMaps,2),allenMask,'split');
taskMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'motor'}),strcmpi(modTypes,'shOther')),2),allenMask,'split');
motorMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'motor'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
% cRange = prctile(fullMap(:),99);
cRange = 0.3;

figure
% subplot(1,3,1)
% mapImg = imshow(fullMap,[0 cRange]);axis image;
% set(mapImg,'AlphaData',~isnan(fullMap)); %make NaNs transparent.
% title('Full model- R^2')
% colormap(mapImg.Parent,jet(256));

subplot(1,2,1)
mapImg = imshow(taskMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(taskMap)); %make NaNs transparent.
title('Task model- R^2')
colormap(mapImg.Parent,inferno(256));

subplot(1,2,2)
mapImg = imshow(motorMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(motorMap)); %make NaNs transparent.
title('Motor model- R^2')
colormap(mapImg.Parent,inferno(256));

%% make movie for task- and motor explained variance
% cMovie = arrayShrink(squeeze(nanmean(fullMovies - taskMovies,2)), allenMask, 'split');
% cMovie = arrayShrink(squeeze(nanmean(fullMovies - motorMovies,2)), allenMask, 'split');
% cMovie = arrayShrink(squeeze(nanmean(fullMovies,2)), allenMask, 'split');
% compareMovie(cMovie,'outline',dorsalMaps.edgeOutlineSplit);

%% single regressor effects
handleMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'rGrab'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
stimMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'rVisStim'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
expMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'expReward'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
novMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,strcmpi(altRecLabels,{'novReward'}),strcmpi(modTypes,'shCurrent')),2),allenMask,'split');

figure
subplot(1,4,1)
cRange = prctile(handleMap(:),99);
mapImg = imshow(handleMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(handleMap)); %make NaNs transparent.
title('Right handle- R^2')
colormap(mapImg.Parent,inferno(256));

subplot(1,4,2)
cRange = prctile(handleMap(:),99);
mapImg = imshow(stimMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(stimMap)); %make NaNs transparent.
title('Right stim - R^2')
colormap(mapImg.Parent,inferno(256));

subplot(1,4,3)
cRange = prctile(expMap(:),99);
mapImg = imshow(expMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(expMap)); %make NaNs transparent.
title('Expert reward - R^2')
colormap(mapImg.Parent,inferno(256));

subplot(1,4,4)
cRange = prctile(expMap(:),99);
mapImg = imshow(novMap,[0 cRange]);axis image;
set(mapImg,'AlphaData',~isnan(novMap)); %make NaNs transparent.
title('Novice reward- R^2')
colormap(mapImg.Parent,inferno(256));


%% R2 reduction for all regressors
clear cRegMod cRegRedM cData
cogLabels = [cogLabels {'expReward'} {'novReward'}];
fullM = aRsq(ismember(fullRecLabels,{'full'}), :, 1); %full model
% cInd = dangerMice;
cInd = true(1,length(animals));

cData  = NaN(2,length(fullRecLabels),sum(cInd));
for x = 1:length(fullRecLabels)
    
    cRegMod = aRsq(x, cInd, ismember(modTypes,'shOther')); %this is the current regressor model. All regressor information.
    cRegRedM = fullM(cInd) - aRsq(x, cInd, ismember(modTypes,'shCurrent')); %this is the current regressor. This has non-redundant information.
    
    cData(1,x,:) = cRegMod;
    cData(2,x,:) = cRegRedM;
    
end

idx = zeros(1,length(fullRecLabels));
idx(ismember(fullRecLabels,otherMotorLabels(~ismember(otherMotorLabels,{'face', 'spontMotor' ,'opMotor'})))) = 1;
idx(ismember(fullRecLabels,cogLabels)) = 2;
idx(ismember(fullRecLabels,sensorLabels)) = 3;
% idx(isnan(cData(1,:,1))) = 0;

figure
[ax, idxGroup] = regressorBoxPlot(squeeze(cData(1,idx>0,:))',fullRecLabels(idx>0),5,subplot(2,2,1:2),[0 1 0],idx(idx>0),0.6);
ax = regressorBoxPlot(squeeze(cData(2,idx>0,:))'.*3,fullRecLabels(idx>0),5,ax,[25 111 61]/255,idxGroup,0.3);
legend(ax.Children([3 1]),{'All information' 'Non-redundant information'},'location','northwest')
title('Regressor comparison');
ylabel('crossVal R^2');
xlim([0 sum(idx>0)+1]);

%% overview figure for all unique R2 maps
rightHs = find(ismember(dorsalMaps.sidesSplit,'R')); %index for labels on the left HS
idx = zeros(1,length(altRecLabels));
idx(ismember(altRecLabels,cogLabels)) = 1;
idx(ismember(altRecLabels,sensorLabels)) = 1;
idx(ismember(altRecLabels,otherMotorLabels(~ismember(otherMotorLabels,{'face', 'spontMotor' ,'opMotor'})))) = 2;

for iMods = 1 : 2
    figure('name',['iMods = ' int2str(iMods)])
    cIdx = find(idx == iMods);
    Cnt = 0;
    
    for x = cIdx
        
        Cnt = Cnt +1;
        cMap = fullMap - arrayShrink(nanmean(corrMaps(:,:,x,strcmpi(modTypes,'shCurrent')),2),allenMask,'split');
        allScale{iMods}(Cnt) = max(abs(prctile(cMap(:),[1 99])));
        
        subplot(ceil(length(cIdx) / 4), 4 ,Cnt);
        
        cRange = max(abs(prctile(cMap(:),[1 99])));
        mapImg = imshow(cMap,[0 cRange]);
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        colormap(mapImg.Parent,inferno(256));hold on;
        title(altRecLabels{x}); colorbar
        drawnow;
        
        hold(mapImg.Parent, 'on');
        for x = 1: length(rightHs)
            plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
        end
        
    end
end

%% R2 reduction for regressor groups
cInd = true(1,animalCnt);

fullM = aRsq(ismember(fullRecLabels,{'full'}), cInd, 1); %full model

spMotorM = aRsq(ismember(fullRecLabels,{'spontMotor'}), cInd, ismember(modTypes,'shOther')); %this is the spont. motor only model. Shuffle everything but motor regressors
spMotorRedM = fullM - aRsq(ismember(fullRecLabels,{'spontMotor'}), cInd, ismember(modTypes,'shCurrent')); %this is the non-redundant spont. motor only model. Shuffle all motor regressors and substract from full.

opMotorM = aRsq(ismember(fullRecLabels,{'opMotor'}), cInd, ismember(modTypes,'shOther')); %this is the operant motor only model. Shuffle everything but motor regressors
opMotorRedM = fullM - aRsq(ismember(fullRecLabels,{'opMotor'}), cInd, ismember(modTypes,'shCurrent')); %this is the non-redundant operant motor only model. Shuffle all motor regressors and substract from full.

motorM = aRsq(ismember(fullRecLabels,{'motor'}), cInd, ismember(modTypes,'shOther')); %this is the motor only model. Shuffle everything but motor regressors
motorRedM = fullM - aRsq(ismember(fullRecLabels,{'motor'}), cInd, ismember(modTypes,'shCurrent')); %this is the non-redundant motor only model. Shuffle all motor regressors and substract from full.

taskM = aRsq(ismember(fullRecLabels,{'motor'}), cInd, ismember(modTypes,'shCurrent')); %this is the task only model. Shuffle only motor regressors to leave task regressors
taskRedM = fullM - motorM; %this is the non-redundant task only model.

clear cData cError
cData(1,:,:) = [taskM' motorM' spMotorM' opMotorM' fullM']; %all reg information
cData(2,:,:) = [taskRedM' motorRedM' spMotorRedM' opMotorRedM' fullM']; %all reg information

figure
[ax, idxGroup] = regressorBoxPlot(squeeze(cData(1,:,:)), {'Task' 'Movement' 'spontMove' 'opMove' 'Full'}, 5, subplot(1,1,1), [0 1 0], [],0.6);
ax = regressorBoxPlot(squeeze(cData(2,:,:)), {'Task' 'Movement' 'spontMove' 'opMove' 'Full'}, 5, ax, [25 111 61]/255, idxGroup,0.6);
legend(ax.Children([3 1]),{'All information' 'Non-redundant information'},'location','northwest')
title('Regressor groups comparison');
ylabel('crossVal R^2'); axis square
xlim([0 6]); axis square
hline(mean(fullM),'k--')
ylim([0 0.5]);

%% Task-dependent R2 reduction for motor regressors
taskM = aRsq(ismember(fullRecLabels,{'motor'}), :, ismember(modTypes,'shCurrent')); %this is the task only model. Shuffle only motor regressors to leave task regressors

cData  = NaN(2,length(otherMotorLabels)+1,animalCnt);
for x = 1:length(otherMotorLabels)
    
    cRegMod = aRsq(ismember(fullRecLabels,otherMotorLabels{x}), :, ismember(modTypes,'shOther')); %this is the current regressor model. All regressor information.
    cRegRedM = aRsq(ismember(fullRecLabels,otherMotorLabels{x}), :, ismember(modTypes,'shOtherMotor')) - taskM; %this is the current regressor + task model. This has non-task specific regressor information.
    
    cData(1,x,:) = cRegMod;
    cData(2,x,:) = cRegRedM;
    
end

idx = double(~ismember(otherMotorLabels,{'face', 'spontMotor' ,'opMotor'}));
figure
[ax, allInd] = regressorPlot(squeeze(cData(1,idx>0,:))',otherMotorLabels(idx>0),5,subplot(2,2,1:2),[0,191,255]/255,idx(idx>0),0.6);
ax = regressorPlot(squeeze(cData(2,idx>0,:))',otherMotorLabels(idx>0),5,ax,[0,0,139]/255,allInd,0.6);
legend(ax.Children([3 1]),{'Task-dependent information' 'Task-independent information'},'location','northwest')
title('Regressor comparison - Task-dependancy');
ylabel('crossVal R^2');
xlim([0 sum(idx>0)+1]);

end

