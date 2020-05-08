function [recV, dType, dataPath, allOpts, allModIdx, allSideIdx, alignIdx, baseLength, frames, stimTimes] = BpodImager_motorReconstruct(cPath,cMod,dType,combOnly,fileExt,getPartModel)
% code to create reconstructed imaging data based on motor regressors. cMod
% informs over which animals to use ('Visual', 'Audio' or 'All'), dType is a
% cell vector containing the labels of the motor regressors to be used for
% reconstruction.
% recV is a cell vector length dType + 1. It contains an averaged
% reconstructed V for each motor regressor. The last cell contains
% reconstruction for all regressors combined. Each cell contains a 4D matrix
% of size frames x dimensions x modality x animals.
% The results need to be convoled with the correct U later to get full
% maps back.
% getPartModel: flag to indicate that averages should come from
% reconstructions that are based on a task,opMotor or spMotor - only
% model.

if ~exist('fileExt','var') || isempty(fileExt)
    fileExt = ''; %this is to be able to use other kinds of models. Default is zero-mean model data.
end

if ~exist('getPartModel','var') || isempty(getPartModel)
    getPartModel = false; %if true, this will use 'Vc' from a a partial-model reconstruction instead of raw data
end

%% select data sets
[dataOverview, motorLabels] = delayDecRecordings;
if ~exist('combOnly','var') || isempty(combOnly)
    combOnly = false;  % flag to only produce data reconstruction from combined regressors. 
end

if getPartModel
   combOnly = true; %This has to be true when asking for a partial-model reconstruction
end

%% general variables
Paradigm = 'SpatialDisc';
% cPath = 'H:\BpodImager\Animals\'; %Widefield data path
% cPath = 'U:\space_managed_data\BpodImager\Animals\'; %Widefield data path on grid server
        
%%
animals = dataOverview(:,1);
Cnt = 0;
baseLength = inf;
postLength = inf;
taskOnly = false;
opMotorOnly = false;
spMotorOnly = false;


for iAnimals = 1:length(animals)
% for iAnimals = 7
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
            if strcmpi(dType,'all')
                dType = recLabels;
            elseif strcmpi(dType,'motor')
                dType = motorLabels;
            elseif strcmpi(dType,'task')
                dType = recLabels(~ismember(recLabels, motorLabels));
                taskOnly = true;
            elseif strcmpi(dType,'opMotor')
                dType = {'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick'};
                opMotorOnly = true;
            elseif strcmpi(dType,'spMotor')
                dType = motorLabels(~ismember(motorLabels,{'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick'}));
                spMotorOnly = true;
            end
            
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
        
        if getPartModel && taskOnly
            load([fPath 'interpVtask.mat'], 'Vtask'); %load model reconstruction data
            Vm = Vtask;
        end
        if getPartModel && spMotorOnly
            load([fPath 'interpVspontMotor.mat'], 'VspontMotor'); %load model reconstruction data
            Vm = VspontMotor;
        end
        
        if getPartModel && opMotorOnly
            load([fPath 'interpVopMotor.mat'], 'VopMotor'); %load model reconstruction data
            Vm = VopMotor;
        end
                
        %% get frametimes to determine #frames / trial and construct modality indices
        load([fPath 'Vc.mat'],'bTrials');
        load([fPath 'interpVc.mat'],'frames');

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
        for iRegs = 1:regRuns
            if getPartModel %use exisiting model reconstruction
                Vm = bsxfun(@minus, Vm, mean(Vm,2));
                Vm = Vm(:,alignIdx{Cnt});
            else
                cInd = ismember(recIdx(~idx), find(ismember(recLabels,dType{iRegs}))); %find current regressors
                data = fullR(:, cInd); %current regressors
                cBeta = dimBeta(cInd,:); %current weights
                Vm = (data * cBeta)'; %model Vc data
                Vm = bsxfun(@minus, Vm, mean(Vm,2));
            end

            if Cnt == 1
                if strcmpi(cMod,'all')
                    recV{iRegs} = NaN(size(Vm,1), frames, 3, length(animals));
                else
                    recV{iRegs} = NaN(size(Vm,1), frames, 3, sum(ismember(dataOverview(:,2), cMod)));
                end
            end
            for iMod = 1:3 %different modalities, 1 is corr. vision, 2 is corr. audio, 3 is all corr. trials. 4:6 is the same thing but for all trials.
                temp = reshape(Vm(:,modIdx(iMod,:)),size(Vm,1),frames,[]);
                recV{iRegs}(:,:,iMod,Cnt) = mean(temp,3); %reconstructed V for current motor regressor and modality
            end
%             temp = reshape(Vm(:,modIdx(6,:)),size(Vm,1),frames,[]); %get trial-by-trial average for all trials
%             trialAvg{Cnt}(iRegs,:,:) = squeeze(mean(temp,2));
        end

        for iMod = 1:6
            allModIdx{Cnt,iMod} = modIdx(iMod,:);
            allSideIdx{Cnt,iMod} = sideIdx(iMod,:);
        end
        clear sideIdx modIdx Vm
    end
end
end