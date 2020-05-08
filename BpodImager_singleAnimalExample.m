%% basic variables
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
trialSegments = [0 55 80 130 160]; %segments to create maps for whole trial regressors
baseLength = inf;
postLength = inf;
motorIdx = 16; %index for zero-lag motor regressor
baseRange = 1:15; %baseline frames
areaIdx = {'MOs'}; %M2
load allenDorsalMapSM
allenMask = dorsalMaps.allenMask;
rightHs = ismember(dorsalMaps.sidesSplit,'R'); %index for labels on the left HS

%% get session data
[dataOverview, motorLabels, sensorLabels, cogLabels, segIdx, segLabels, segIdxRealign] = delayDecRecordings;
sensorLabels = [sensorLabels cogLabels]; %combine sensory and cognitive labels, so sensory is really 'non-motor' or 'task' here
segIdxRealign{2} = 46:65;

recs = dataOverview(:,3)';
animals = dataOverview(:,1)';
targRec = 7; %target recording (in the paper it is 7 / 22 recordings)

fPath = [cPath animals{targRec} filesep 'SpatialDisc' filesep recs{targRec} filesep]; %Widefield data path
load([fPath 'regData.mat'],'fullR','recIdx','idx','recLabels', 'trialIdx');
load([fPath 'opts2.mat'], 'opts'); %load opts
load([fPath 'interpVc.mat'],'frames');
load([fPath 'Vc.mat'],'U','bTrials');

U = alignAllenTransIm(U,opts.transParams);
U = arrayShrink(U(:,1:size(allenMask,2),:),allenMask,'merge');
U = arrayShrink(U,allenMask,'split');

% load behavior
bhvFile = strsplit(fPath,filesep);
bhvFile = dir([fPath bhvFile{end-3} '*' bhvFile{end-2} '*.mat']);
load([fPath bhvFile.name]);
        
%% realign data so baseline is aligned to handle grab and poststim to stimulus
trialIdx = unique(ceil(find(~trialIdx)/frames)); %rebuild as single trial index (instead of frames)
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
alignIdx = sort([baseIdx;stimIdx]);

% align Vc
frames = baseLength + postLength;
Vc = Vc(:,alignIdx);

% get trial indicices
sucInd = SessionData.Rewarded(bTrials) & SessionData.Assisted(bTrials); %find succesful unisensory trials
modIdx(1,:) = reshape(repmat(SessionData.StimType(bTrials) == 1 & sucInd,frames,1),[],1); % correct visual trials
modIdx(2,:) = reshape(repmat(SessionData.StimType(bTrials) == 2 & sucInd,frames,1),[],1); % correct audio trials
modIdx(3,:) = sum(modIdx(1:2,:)); % all correct trials
modIdx = reshape(modIdx,size(modIdx,1),frames,[]);
modIdx = modIdx(:,:,trialIdx); %remove non-used trials
modIdx = reshape(modIdx,size(modIdx,1),[]);

% get area index and compute traces
ind = ismember(dorsalMaps.labelsSplit,areaIdx) & rightHs;
areaCoord = ~poly2mask(dorsalMaps.edgeOutlineSplit{ind}(:,2), dorsalMaps.edgeOutlineSplit{ind}(:,1),size(allenMask,1),size(allenMask,2));
areaU = arrayShrink(U, areaCoord, 'merge'); %piexels in target area
cData = mean(areaU(~isnan(areaU(:,1)),:),1) * Vc;
cData = reshape(cData,frames,[]);
