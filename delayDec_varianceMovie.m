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

%%
Cnt = 0;
baseLength = inf;
postLength = inf;
taskMovie = NaN(sum(~allenMask(:)), 179, length(animals), 'single');
motorMovie = NaN(sum(~allenMask(:)), 179, length(animals), 'single');
fullMovie = NaN(sum(~allenMask(:)), 179, length(animals), 'single');

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
        load([fPath 'interpVc.mat'],'Vc','frames');
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
        Vc = Vc(:,alignIdx{Cnt}); %reduce fullR to only include aligned baseline and poststim data
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
        
        %% cycle through frames and reconstruct each one
        Vc = reshape(Vc,size(Vc,1),[]);
        U = U(:,1:size(Vc,1));
        
        cInd = ismember(recIdx(~idx), find(~ismember(recLabels,motorLabels))); %find task regressors
        cInd(std(fullR,0,1) == 0) = false; %double check that there are no empty regressors after realignment
        [~, cBeta] = ridgeMML(Vc', fullR(:, cInd), true); %get ridge penalties and beta weights.
        Vt = (fullR(:, cInd) * cBeta)'; %model task data
        Vt = reshape(Vt,size(Vt,1),frames,[]);

        cInd = ismember(recIdx(~idx), find(ismember(recLabels,motorLabels))); %find motor regressors
        cInd(std(fullR,0,1) == 0) = false; %double check that there are no empty regressors after realignment
        [~, cBeta] = ridgeMML(Vc', fullR(:, cInd), true); %get ridge penalties and beta weights.
        Vm = (fullR(:, cInd) * cBeta)'; %model movement data                
        Vm = reshape(Vm,size(Vm,1),frames,[]);
        Vc = reshape(Vc,size(Vc,1),frames,[]);
        
        Vf = (fullR * dimBeta)'; %full model       
        Vf = reshape(Vf,size(Vf,1),frames,[]);  
        
        for iFrames = 1 : frames
            
            cFrame = squeeze(Vc(:,iFrames,:));
            mFrame = squeeze(Vm(:,iFrames,:));
            tFrame = squeeze(Vt(:,iFrames,:));
            fFrame = squeeze(Vf(:,iFrames,:));
            
            covVc = cov(cFrame');  % S x S
            covVm = cov(mFrame');  % S x S
            covVt = cov(tFrame');  % S x S
            covVf = cov(fFrame');  % S x S
        
            mCovV = bsxfun(@minus, mFrame, mean(mFrame,2)) * cFrame' / (size(cFrame, 2) - 1);  % S x S
            tCovV = bsxfun(@minus, tFrame, mean(tFrame,2)) * cFrame' / (size(cFrame, 2) - 1);  % S x S
            fCovV = bsxfun(@minus, fFrame, mean(fFrame,2)) * cFrame' / (size(cFrame, 2) - 1);  % S x S

            covPm = sum((U * mCovV) .* U, 2)';  % 1 x P
            covPt = sum((U * tCovV) .* U, 2)';  % 1 x P
            covPf = sum((U * fCovV) .* U, 2)';  % 1 x P
            
            varVc = sum((U * covVc) .* U, 2)';  % 1 x P
            varVm = sum((U * covVm) .* U, 2)';  % 1 x P
            varVt = sum((U * covVt) .* U, 2)';  % 1 x P
            varVf = sum((U * covVf) .* U, 2)';  % 1 x P
            
            stdVm = varVc .^ 0.5 .* varVm .^ 0.5; % 1 x P
            stdVt = varVc .^ 0.5 .* varVt .^ 0.5; % 1 x P
            stdVf = varVc .^ 0.5 .* varVf .^ 0.5; % 1 x P

            mCorrMap = (covPm ./ stdVm).^2';
            tCorrMap = (covPt ./ stdVt).^2';
            fCorrMap = (covPf ./ stdVf).^2';

            motorMovie(:,iFrames,Cnt) = fCorrMap - tCorrMap; %unique motor variance (full - task)
            taskMovie(:,iFrames,Cnt) = fCorrMap - mCorrMap; % unique task variance (full - motor)
            fullMovie(:,iFrames,Cnt) = fCorrMap; % all explained variance
        end
        clear sideIdx modIdx Vm
        
        if rem(iAnimals, round(length(animals)/5)) == 0
            fprintf(1, 'Done. Current recording is %d out of %d\n', iAnimals, length(animals));
            toc;
        end
    end
end

%% different R2 movies
compareMovie(arrayShrink(nanmean(taskMovie,3),allenMask,'split'));
compareMovie(arrayShrink(nanmean(motorMovie,3),allenMask,'split'));
compareMovie(arrayShrink(nanmean(fullMovie,3),allenMask,'split'));
