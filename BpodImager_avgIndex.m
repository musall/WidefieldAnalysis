function [dataPath, allOpts, allModIdx, allSideIdx, alignIdx, baseLength, frames, stimTimes] = BpodImager_avgIndex(cMod, combOnly, fileExt)
% code to get different indices to create PSTHs from imaging data.
% dataPath is path to each dataset. allOpts is for alignment to the CCF.
% allModIdx is a cell that contains indices for either visual, auditory, or
% both kind of trials for each animal. 
% allSideIdx is the as allModIdx for left and right. 
% alignIdx can be used to the subset of continious data that yields handle
% and stimulus aligned PSTHs. 
% baseLength is the duration of the new baseline in seconds (relative to
% first possible stimulus onset in all trials).
% frames is the number of frames in double-aligned PSTHs.
% stimTimes containes the onset of all stimuli, relative to trialOnset in
% seconds.

if ~exist('fileExt','var') || isempty(fileExt)
    fileExt = ''; %this is to be able to use other kinds of models. Default is zero-mean model data.
end

%% select data sets
dataOverview = delayDecRecordings;
if ~exist('combOnly','var') || isempty(combOnly)
    combOnly = false;  % flag to only produce data reconstruction from combined regressors
end

%% general variables
Paradigm = 'SpatialDisc';
cPath = 'U:\space_managed_data\BpodImager\Animals\'; %Widefield data path on grid server
sPath = 'U:\space_managed_data\BpodImager\Animals\'; %Widefield data path on grid server
        
%%
animals = dataOverview(:,1);
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
        load([fPath fileExt 'regData.mat'],'fullR','trialIdx','recIdx','idx','recLabels'); %load model data
        allOpts(Cnt) = load([fPath 'opts2.mat']); %load opts

        %% get frametimes to determine #frames / trial and construct modality indices
        load([fPath 'Vc.mat'],'bTrials');
        load([fPath 'interpVc.mat'],'frames');
        trialIdx = unique(ceil(find(~trialIdx)/frames)); %rebuild as single trial index (instead of frames)

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
        
        %% return indices
        for iMod = 1:6
            allModIdx{Cnt,iMod} = modIdx(iMod,:);
            allSideIdx{Cnt,iMod} = sideIdx(iMod,:);
        end
        clear sideIdx modIdx Vm
    end
end
end
 
        