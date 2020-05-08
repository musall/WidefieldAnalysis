%code to produce dPrime maps for novice versus expert decision when
%deconstructing data into motor and task.

cMod = 'all';
if ~exist('fileExt','var') || isempty(fileExt)
    fileExt = ''; %this is to be able to use other kinds of models. Default is zero-mean model data.
end

[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, segIdxRealign] = delayDecRecordings;

if ~exist('combOnly','var') || isempty(combOnly)
    combOnly = false;  % flag to only produce data reconstruction from combined regressors
end
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts
animals = dataOverview(:,1);
tic

%% general variables
Paradigm = 'SpatialDisc';
% cPath = 'H:\BpodImager\Animals\'; %Widefield data path
cPath = 'U:\space_managed_data\BpodImager\Animals\'; %Widefield data path on grid server

% %get allen maps
load('allenDorsalMap.mat')
addEdgeOutlinesToDorsalMap(dorsalMaps); %make area figure for orientation
mask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
[x1, y1] = size(mask);
fPath = [cPath animals{1} filesep Paradigm filesep dataOverview{1,3} filesep];
load([fPath 'snapshot_1.mat'])
[x2, y2] = size(snap);
allenMask = mask(1:min([x1 x2]), 1:min([y1 y2])); %cut mask to size
rightHs = find(ismember(dorsalMaps.sidesSplit,'R')); %index for labels on the left HS

cSegs = cat(2,segIdxRealign{ismember(segLabels,{'Stim1' 'Stim2'})}); %time index for stimulus segments
cSegs = [cSegs(1) segIdxRealign{ismember(segLabels,{'Wait'})}(round(end/2))]; %add delay period
cSegs = cSegs([1 end]);

%%
Cnt = 0;
baseLength = inf;
postLength = inf;

for iAnimals = 1:length(animals)
    if strcmpi(dataOverview{iAnimals,2},cMod) || strcmpi(cMod,'all')        
        %% load data
        Cnt = Cnt +1;
        fPath = [cPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
        dataPath{Cnt} = fPath; %store current data path so it is easier to load in U later
        
        load([fPath 'Snapshot_1.mat']);
        load([fPath 'mask.mat'])
        load([fPath fileExt 'dimBeta.mat'],'dimBeta');
        load([fPath 'regData.mat'],'fullR','trialIdx','recIdx','idx','recLabels'); %load model data
        allOpts(Cnt) = load([fPath 'opts2.mat']); %load opts
        
        %% check dType and rebuild if required
        if Cnt == 1
            
            dType = recLabels;
            if ischar(dType)
                dType = cellstr(dType);
            end
            
            if combOnly  %if only combined results are requested, combine regressors and only run once
                if iscellstr(dType)
                    dType{1} = dType;
                else
                    dType{1} = cat(2,dType{:});
                end
                dType(2:end) = [];
                regRuns = 1;
            else
                regRuns = length(dType);
            end
        end
        
        %% get frametimes to determine #frames / trial and construct modality indices
        load([fPath 'interpVc.mat'],'frames');
        load([fPath 'Vc.mat'],'bTrials', 'U');
        U = alignAllenTransIm(U,allOpts(iAnimals).opts.transParams);
        U = arrayShrink(U(1:size(allenMask,1),1:size(allenMask,2),:), allenMask);

        % load behavior and get modality indices
        bhvFile = strsplit(fPath,filesep);
        bhvFile = dir([fPath bhvFile{end-3} '*' bhvFile{end-2} '*.mat']);
        load([fPath bhvFile.name]);
        
        for iTrials = 1 : length(SessionData.Rewarded)
            try
                leverTimes = [reshape(SessionData.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
                    reshape(SessionData.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
                    reshape(SessionData.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
                stimGrab = leverTimes(find(leverTimes == SessionData.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
                stimTimes{Cnt}(iTrials) = SessionData.RawEvents.Trial{iTrials}.Events.Wire3High - stimGrab; %time of stimulus onset - measured from soundcard
            catch
                stimTimes{Cnt}(iTrials) = NaN;
            end
        end
        clear stimGrab leverTimes
        
        %% realign data so baseline is aligned to handle grab and poststim to stimulus
        stimOn = sum(fullR(:,ismember(recIdx(~idx),find(ismember(recLabels,{'lVisStim' 'rVisStim' 'lAudStim' 'rAudStim'})))),2); %index for stimulus onset in all trials
        stimOn = find([0;diff(stimOn)] > 0.5) - 1;

        % index for baseline (time before first possible stimulus onset)
        baseLength = min([min(unique(rem(stimOn,frames)))-1 baseLength]);
        baseIdx = repmat((0:frames:size(fullR,1)-1)',1,baseLength);
        baseIdx = bsxfun(@plus,baseIdx, 1 : baseLength);
        baseIdx = baseIdx(:);

        % index for post stimulus time
        postLength = min([frames - max(unique(rem(stimOn,frames))) postLength]); %shortest possible poststim duration
        stimIdx = repmat(stimOn, 1, postLength);
        stimIdx = bsxfun(@plus,stimIdx, 0 : postLength-1);
        stimIdx = stimIdx(:);
        alignIdx{Cnt} = sort([baseIdx;stimIdx]);
        
        if ~isempty(fileExt)
            load([fPath fileExt 'regData.mat'],'fullR','trialIdx','recIdx','idx','recLabels'); %load model data
        end
        
        trialIdx = unique(ceil(find(~trialIdx)/frames)); %rebuild as single trial index (instead of frames)
        fullR = bsxfun(@minus, fullR, mean(fullR, 1)); %make sure design matrix is zero-mean
        fullR = fullR(alignIdx{Cnt},:); %reduce fullR to only include aligned baseline and poststim data
        frames = postLength + baseLength; %new single trial duration in frames
        
        %% build lots of indices
        sucInd = SessionData.Rewarded(bTrials) & SessionData.Assisted(bTrials); %find succesful unisensory trials
        modIdx(1,:) = reshape(repmat(SessionData.StimType(bTrials) == 1 & sucInd,frames,1),[],1); % correct visual trials
        modIdx(2,:) = reshape(repmat(SessionData.StimType(bTrials) == 2 & sucInd,frames,1),[],1); % correct audio trials
        modIdx(3,:) = sum(modIdx(1:2,:)); % all correct trials
        modIdx(4,:) = reshape(repmat(SessionData.StimType(bTrials) == 1, frames,1),[],1); % all visual trials
        modIdx(5,:) = reshape(repmat(SessionData.StimType(bTrials) == 2, frames,1),[],1); % all audio trials        
        modIdx(6,:) = sum(modIdx(4:5,:)); % all trials
        modIdx = reshape(modIdx,size(modIdx,1),frames,[]);
        modIdx = modIdx(:,:,trialIdx); %remove non-used trials
        modIdx = reshape(modIdx,size(modIdx,1),[]);

        sideIdx(1,:) = reshape(repmat(SessionData.CorrectSide(bTrials) == 1 & sucInd,frames,1),[],1); % correct left trials
        sideIdx(2,:) = reshape(repmat(SessionData.CorrectSide(bTrials) == 2 & sucInd,frames,1),[],1); % correct right trials
        sideIdx(3,:) = sum(sideIdx(1:2,:)); % all correct trials
        sideIdx(4,:) = reshape(repmat(SessionData.CorrectSide(bTrials) == 1, frames,1),[],1); % correct left trials
        sideIdx(5,:) = reshape(repmat(SessionData.CorrectSide(bTrials) == 2, frames,1),[],1); % correct right trials
        sideIdx(6,:) = sum(sideIdx(4:5,:)); % all trials    
        sideIdx = reshape(sideIdx,size(sideIdx,1),frames,[]);
        sideIdx = sideIdx(:,:,trialIdx); %remove non-used trials
        sideIdx = reshape(sideIdx,size(sideIdx,1),[]);
        
        %% cycle through regressors and reconstruct each one
        x = [3 5]; % segments from segIdx for d' calculation
        for iRegs = 1:2 %1 is task, 2 is movement
            
            if iRegs == 1
                cInd = ismember(recIdx(~idx), find(~ismember(recLabels,motorLabels))); %find task regressors
            elseif iRegs == 2
                cInd = ismember(recIdx(~idx), find(ismember(recLabels,motorLabels))); %find motor regressors
            end
            
            data = fullR(:, cInd); %current regressors
            cBeta = dimBeta(cInd,:); %current weights
            Vm = (data * cBeta)'; %model Vc data
            Vm = reshape(Vm,size(Vm,1),frames,[]);
            U = U(:,1:size(cBeta,2));
                
            for iSegs = 1:2 %stim1 and delay window
                
                % compute standard deviation for expert and novice decisions (all correct trials)
                cData = reshape(Vm(:,modIdx(3,:)),size(Vm,1),frames,[]);
                cData = squeeze(mean(cData(:,segIdxRealign{x(iSegs)},:),2));
                covV = cov(cData');  % S x S
                varP = sum((U * covV) .* U, 2);  % 1 x P
                stdP = sqrt(varP); %standard deviation map for correct trials
                
                % compute average for expert and novice decisions
                visData = reshape(Vm(:,modIdx(1,:)),size(Vm,1),frames,[]);
                visData = squeeze(mean(visData(:,segIdxRealign{x(iSegs)},:),2));
                audData = reshape(Vm(:,modIdx(2,:)),size(Vm,1),frames,[]);
                audData = squeeze(mean(audData(:,segIdxRealign{x(iSegs)},:),2));
                
                if audExp(iAnimals)
                    expD = mean(U * audData,2) - mean(U * visData,2); % audio - vision PSTH
                elseif visExp(iAnimals)
                    expD = mean(U * visData,2) - mean(U * audData,2); % vision - audio PSTH
                end
                expD = expD ./ stdP;
                aExpD{iRegs,iSegs}(:,:,Cnt) = arrayShrink(expD,allenMask,'split');
                
            end
        end
        clear sideIdx modIdx Vm
        
        if rem(iAnimals, round(length(animals)/5)) == 0
            fprintf(1, 'Done. Current recording is %d out of %d\n', iAnimals, length(animals));
            toc;
        end
    end
end

%%
figure('name','task')
subplot(2,1,1)
cMap = nanmean(aExpD{1,1},3);
cRange = max(abs(prctile(cMap(:),[1 99])));
mapImg = imshow(cMap,[-2 2]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
title('Stim1')

subplot(2,1,2)
cMap = nanmean(aExpD{1,2},3);
cRange = max(abs(prctile(cMap(:),[1 99])));
mapImg = imshow(cMap,[-2 2]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
title('Delay')

figure('name','movement')
subplot(2,1,1)
cMap = nanmean(aExpD{2,1},3);
cRange = max(abs(prctile(cMap(:),[1 99])));
mapImg = imshow(cMap,[-1 1]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
title('Stim1')

subplot(2,1,2)
cMap = nanmean(aExpD{2,2},3);
cRange = max(abs(prctile(cMap(:),[1 99])));
mapImg = imshow(cMap,[-1 1]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
title('Delay')





        

 
        