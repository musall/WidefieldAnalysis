function BpodImager_delayRegressModel(cPath,Animal,Rec)

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

%% general variables
Paradigm = 'SpatialDisc';
cPath = [cPath Animal filesep Paradigm filesep Rec filesep]; %Widefield data path
sPath = ['U:\space_managed_data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
% sPath = ['/sonas-hs/churchland/hpc/home/space_managed_data/BpodImager/Animals/' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
sRate = 30; % Sampling rate of imaging in Hz
preStimDur = 1.8;       % Duration of trial before lever grab in seconds
postStimDur = 4.5;      % Duration of trial after lever grab onset in seconds
frames = (preStimDur + postStimDur) * sRate; %nr of frames per trial
trialDur = (frames * (1/sRate)); %duration of trial in seconds

%other variables
mPreTime = 0.5;     % precede motor events to capture preparatory activity in seconds
mPostTime = 1;      % follow motor events for mPostStim in seconds
motorIdx = [-((mPreTime * sRate): -1 : 1) 0 (1:(mPostTime * sRate))]; %index for design matrix to cover pre- and post motor action
tapDur = 0.25;      % minimum time of lever contact, required to count as a proper grab.
piezoLine = 2;      % channel in the analog data that contains data from piezo sensor
stimLine = 6;       % channel in the analog data that contains stimulus trigger.
leverMoveDur = 0.25; %duration of lever movement. this is used to orthogonalize video against lever movement.
leverMoveDur = ceil(leverMoveDur * sRate); %convert to frames

bhvDimCnt = 200;    % number of dimensions from behavioral videos that are used as regressors.
dims = 200; %number of dimensions from V that will be used in the model
gaussShift = 1;     % inter-frame interval between regressors. Will use only every 'gaussShift' regressor and convolve with gaussian of according FHWM to reduce total number of used regressors.

%% load data
if exist([cPath 'Vc.mat'],'file') ~= 2 %check if data file exists and get from server otherwise
    copyfile([sPath 'Vc.mat'],[cPath 'Vc.mat']);
    copyfile([sPath 'mask.mat'],[cPath 'mask.mat']);
    bhvFile = dir([sPath filesep Animal '_' Paradigm '*.mat']);
    copyfile([sPath bhvFile.name],[cPath bhvFile.name]);
end

load([cPath 'mask.mat'])
load([cPath 'Vc.mat'])
Vc = Vc(1:dims,:,:);
U = U(:,:,1:dims);

bhvFile = dir([cPath filesep Animal '_' Paradigm '*.mat']);
load([cPath bhvFile(1).name],'SessionData'); %load behavior data
SessionData.TrialStartTime = SessionData.TrialStartTime * 86400; %convert trailstart timestamps to seconds

% ensure there are not too many trials in Vc
ind = trials > SessionData.nTrials;
trials(ind) = [];
Vc(:,:,ind) = [];
trialCnt = length(trials);

if ~exist('bTrials','var')
    bTrials = trials;
end
bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dataset

%% load behavior data
if exist([cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],'file') ~= 2 %check if svd behavior exists on hdd and pull from server otherwise
    if exist([cPath 'BehaviorVideo' filesep]) ~= 2
        mkdir([cPath 'BehaviorVideo' filesep]);
    end
    copyfile([sPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],[cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'],[cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'FilteredPupil.mat'],[cPath 'BehaviorVideo' filesep 'FilteredPupil.mat']);
end

load([cPath 'BehaviorVideo' filesep 'SVD_CombinedSegments.mat'],'vidV'); %load behavior video data
V1 = vidV(:,1:bhvDimCnt); %behavioral video regressors
load([cPath 'BehaviorVideo' filesep 'motionSVD_CombinedSegments.mat'],'vidV'); %load abs motion video data
V2 = vidV(:,1:bhvDimCnt); % motion regressors

load([cPath 'BehaviorVideo' filesep 'FilteredPupil.mat']); %load pupil data
%check if timestamps from pupil data are shifted against bhv data
timeCheck1 = (SessionData.TrialStartTime(1)) - (pTime{1}(1)); %time difference between first acquired frame and onset of first trial
timeCheck2 = (SessionData.TrialStartTime(1)) - (bTime{1}(1)); %time difference between first acquired frame and onset of first trial
if (timeCheck1 > 3590 && timeCheck1 < 3610) && (timeCheck2 > 3590 && timeCheck2 < 3610) %timeshift by one hour (+- 10seconds)
    warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
    for iTrials = 1 : length(pTime)
        pTime{iTrials} = pTime{iTrials} + 3600; %add one hour
        bTime{iTrials} = bTime{iTrials} + 3600; %add one hour
    end
elseif timeCheck1 > 30 || timeCheck1 < -30 || timeCheck2 > 30 || timeCheck2 < -30
    error('Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds.')
end

    
%% find events in BPod time - All timestamps are relative to stimulus onset event to synchronize to imaging data later
% pre-allocate vectors
lickL = cell(1,trialCnt);
lickR = cell(1,trialCnt);
leverIn = NaN(1,trialCnt);
levGrabL = cell(1,trialCnt);
levGrabR = cell(1,trialCnt);
levReleaseL = cell(1,trialCnt);
levReleaseR = cell(1,trialCnt);
water = NaN(1,trialCnt);

for iTrials = 1:trialCnt
        
    leverTimes = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
    stimGrab(iTrials) = leverTimes(find(leverTimes == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
    
    try
        stimTime(iTrials) = bhv.RawEvents.Trial{iTrials}.Events.Wire3High - stimGrab(iTrials); %time of stimulus onset - measured from soundcard
    catch
        stimTime(iTrials) = NaN;
    end
    
    %check for spout motion
    if isfield(bhv.RawEvents.Trial{iTrials}.States,'MoveSpout') 
        spoutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1) - stimGrab(iTrials);
        
        %also get time when the other spout was moved out at 
        if bhv.Rewarded(iTrials)
            spoutOutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.Reward(1) - stimGrab(iTrials);
        else
            spoutOutTime(iTrials) = bhv.RawEvents.Trial{iTrials}.States.HardPunish(1) - stimGrab(iTrials);
        end
    else
        spoutTime(iTrials) = NaN;
        spoutOutTime(iTrials) = NaN;
    end
    
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port1In') %check for licks
        lickL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port1In;
        lickL{iTrials}(lickL{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickL{iTrials} = lickL{iTrials} - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port3In') %check for right licks
        lickR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port3In;
        lickR{iTrials}(lickR{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickR{iTrials} = lickR{iTrials} - stimGrab(iTrials);
    end
    
    leverIn(iTrials) = min(bhv.RawEvents.Trial{iTrials}.States.Reset(:)) - stimGrab(iTrials); %first reset state causes lever to move in
        
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2High') %check for left grabs
        levGrabL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2High - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1High') %check for right grabs
        levGrabR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1High - stimGrab(iTrials);
    end
    
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2Low') %check for left release
        levReleaseL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2Low - stimGrab(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1Low') %check for right release
        levReleaseR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1Low - stimGrab(iTrials);
    end
    
    if ~isnan(bhv.RawEvents.Trial{iTrials}.States.Reward(1)) %check for reward state
        water(iTrials) = bhv.RawEvents.Trial{iTrials}.States.Reward(1) - stimGrab(iTrials);
    end
end

maxStimRegs = length(min(round((preStimDur + stimTime) * sRate)) : (preStimDur + postStimDur) * sRate); %maximal number of required stimulus regressors
maxSpoutRegs = length(min(round((preStimDur + spoutTime) * sRate)) : (preStimDur + postStimDur) * sRate); %maximal number of required stimulus regressors

%% build regressors - create design matrix based on event times
%basic time regressors
timeR = logical(diag(ones(1,frames)));

lGrabR = cell(1,trialCnt);
lGrabRelR = cell(1,trialCnt);
rGrabR = cell(1,trialCnt);
rGrabRelR = cell(1,trialCnt);
lLickR = cell(1,trialCnt);
rLickR = cell(1,trialCnt);
leverInR = cell(1,trialCnt);

lVisStimR = cell(1,trialCnt);
rVisStimR = cell(1,trialCnt);
lAudStimR = cell(1,trialCnt);
rAudStimR = cell(1,trialCnt);
spoutR = cell(1,trialCnt);
spoutOutR = cell(1,trialCnt);

visRewardR = cell(1,trialCnt);
audRewardR = cell(1,trialCnt);
prevRewardR = cell(1,trialCnt);

ChoiceR = cell(1,trialCnt);

prevChoiceR = cell(1,trialCnt);
prevModR = cell(1,trialCnt);

waterR = cell(1,trialCnt);
fastPupilR = cell(1,trialCnt);
slowPupilR = cell(1,trialCnt);

whiskR = cell(1,trialCnt);
noseR = cell(1,trialCnt);
piezoR = cell(1,trialCnt);
piezoMoveR = cell(1,trialCnt);
faceR = cell(1,trialCnt);
bodyR = cell(1,trialCnt);

%%
tic
for iTrials = 1:trialCnt
    %% vis/aud stim - regressors cover the remaining trial after stimulus onset
    stimIdx = round((preStimDur + stimTime(iTrials)) * sRate) : (preStimDur + postStimDur)*sRate; %index for which part of the trial should be covered by stim regressors
       
    % vision
    lVisStimR{iTrials} = false(frames, maxStimRegs);
    rVisStimR{iTrials} = false(frames, maxStimRegs);
    if bhv.StimType(iTrials) == 1 || bhv.StimType(iTrials) == 3 %visual or mixed stimulus
        if bhv.CorrectSide(iTrials) == 1
            lVisStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        else
            rVisStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        end
    end
    
    % audio
    lAudStimR{iTrials} = false(frames, maxStimRegs);
    rAudStimR{iTrials} = false(frames, maxStimRegs);
    if bhv.StimType(iTrials) == 2 || bhv.StimType(iTrials) == 3 %auditory or mixed stimulus
        if bhv.CorrectSide(iTrials) == 1
            lAudStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        else
            rAudStimR{iTrials}(:, 1:length(stimIdx)) = timeR(:, stimIdx);
        end
    end
    
    if gaussShift > 1
        % subsample regressors
        lVisStimR{iTrials} = lVisStimR{iTrials}(:,1:gaussShift:end);
        rVisStimR{iTrials} = rVisStimR{iTrials}(:,1:gaussShift:end);
        lAudStimR{iTrials} = lAudStimR{iTrials}(:,1:gaussShift:end);
        rAudStimR{iTrials} = rAudStimR{iTrials}(:,1:gaussShift:end);
    end
    
    %% spout regressors
    spoutIdx = round((preStimDur + spoutTime(iTrials)) * sRate) : (preStimDur + postStimDur)*sRate; %index for which part of the trial should be covered by spout regressors
    spoutR{iTrials} = false(frames, maxSpoutRegs);
    spoutR{iTrials}(:, 1:length(spoutIdx)) = timeR(:, spoutIdx);
    
    spoutOutR{iTrials} = false(frames, 3);
    if ~isnan(spoutOutTime(iTrials))
        spoutOut = round((preStimDur + spoutOutTime(iTrials)) * sRate); %time when opposing spout moved out again
        spoutOutR{iTrials}(spoutOut : spoutOut + 2, :) = diag(ones(1,3));
    end
    
    %% lick regressors
    lLickR{iTrials} = false(frames, length(motorIdx));
    rLickR{iTrials} = false(frames, length(motorIdx));
    
    for iRegs = 0 : length(motorIdx)-1
        licks = lickL{iTrials} - (mPreTime - (iRegs * 1/sRate));
        lLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
        
        licks = lickR{iTrials} - (mPreTime - (iRegs * 1/sRate));
        rLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
    end
    
    if gaussShift > 1
        % subsample regressors
        lLickR{iTrials} = lLickR{iTrials}(:,1:gaussShift:end);
        rLickR{iTrials} = rLickR{iTrials}(:,1:gaussShift:end);
    end      
    
    %% lever in
    leverInR{iTrials} = false(frames, leverMoveDur);
    leverShift = round((preStimDur + leverIn(iTrials))* sRate); %timepoint in frames when lever moved in, relative to lever grab
    
    if ~isnan(leverShift)
        if leverShift > 0 %lever moved in during the recorded trial
            leverInR{iTrials}(leverShift : leverShift + leverMoveDur -1, :) = diag(ones(1, leverMoveDur));
        elseif (leverShift + leverMoveDur) > 0  %lever was moving before data was recorded but still moving at trial onset
            leverInR{iTrials}(1 : leverMoveDur + leverShift, :) = [zeros(leverMoveDur + leverShift, abs(leverShift)) diag(ones(1, leverMoveDur + leverShift))];
        end
    end
    
    if gaussShift > 1
        % subsample regressors
        leverInR{iTrials} = leverInR{iTrials}(:,1:gaussShift:end);
    end

    %% choice and reward
    stimShift = round((stimTime(iTrials))* sRate) - sRate; %timepoint in frames when the stimulus was presented. This is the shift relative to the expectation that the stimulus comes up 1s after grabing the lever.
    
    stimTemp = false(frames,frames);
    if stimShift > 0 %stim came later than 1s from lever grab
        stimTemp(:, 1: frames - stimShift) = timeR(:, stimShift+1:end);
    else %stim came earlier than 1s from lever grab
        stimTemp(:, abs(stimShift) + 1 : end) = timeR(:, 1: frames + stimShift);
    end
    stimTemp(:,end-4:end) = []; %don't use the last timepoint to avoid rank defficient design matrix

    visRewardR{iTrials} = false(size(stimTemp));
    audRewardR{iTrials} = false(size(stimTemp));
    if bhv.Rewarded(iTrials) %rewarded
        if bhv.StimType(iTrials)  == 1 %vision
            visRewardR{iTrials} = stimTemp; %visual trial, rewarded
        elseif bhv.StimType(iTrials)  == 2 %audio
            audRewardR{iTrials} = stimTemp; %audio trial, rewarded
        end
    end
    
    % get L/R choices as binary design matrix
    ChoiceR{iTrials} = false(size(stimTemp));
    if bhv.ResponseSide(iTrials) == 1
        ChoiceR{iTrials} = stimTemp;
    end
          
    % previous trial regressors
    if iTrials == 1 %don't use first trial
        prevRewardR{iTrials} = NaN(size(timeR(:,1:end-4)));
        prevChoiceR{iTrials} = NaN(size(timeR(:,1:end-4)));
        prevModR{iTrials} = NaN(size(timeR(:,1:end-4)));
        
    else %for all subsequent trials
        % same as for regular choice regressors but for prevoious trial
        prevChoiceR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.ResponseSide(bTrials(iTrials)-1) == 1
            prevChoiceR{iTrials} = timeR(:,1:end-4);
        end
             
        prevModR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.StimType(bTrials(iTrials)-1) == 1 || SessionData.StimType(bTrials(iTrials)-1) == 3
            prevModR{iTrials} = timeR(:,1:end-4); % if previous trial was vision
        end
        
        prevRewardR{iTrials} = false(size(timeR(:,1:end-4)));
        if SessionData.Rewarded(bTrials(iTrials)-1) %last trial was rewarded
            prevRewardR{iTrials} = timeR(:,1:end-4);
        end
    end
       
    if gaussShift > 1
        % subsample regressors
        visRewardR{iTrials} = visRewardR{iTrials}(:,1:gaussShift:end);
        audRewardR{iTrials} = audRewardR{iTrials}(:,1:gaussShift:end);
        prevRewardR{iTrials} = prevRewardR{iTrials}(:,1:gaussShift:end);

        ChoiceR{iTrials} = ChoiceR{iTrials}(:,1:gaussShift:end);
        prevChoiceR{iTrials} = prevChoiceR{iTrials}(:,1:gaussShift:end);

        prevModR{iTrials} = prevModR{iTrials}(:,1:gaussShift:end);
    end
    
    %determine timepoint of reward given
    waterR{iTrials} = false(frames, sRate);
    if ~isnan(water(iTrials)) && ~isempty(water(iTrials))
        waterOn = round((preStimDur + water(iTrials)) * sRate); %timepoint in frames when reward was given
        waterR{iTrials}(:, 1: size(timeR,2) - waterOn + 1) = timeR(:, waterOn:end);
    end
    
    if gaussShift > 1
        waterR{iTrials} = waterR{iTrials}(:,1:gaussShift:end); % subsample regressor
    end
    
    %% lever grab / release
    lGrabR{iTrials} = false(frames, length(motorIdx));
    lGrabRelR{iTrials} = false(frames, length(motorIdx));
    
    rGrabR{iTrials} = false(frames, length(motorIdx));
    rGrabRelR{iTrials} = false(frames, length(motorIdx));
    
    [grabL,grabRelL] = checkLevergrab(tapDur,postStimDur,levGrabL{iTrials},levReleaseL{iTrials},1/sRate);
    grabL(grabL < -preStimDur) = []; %ensure there are no grab events too early in the trial
    if ~isempty(grabL);grabRelL(grabRelL<grabL(1)) = []; end %make sure there are no release events before first grab
    
    [grabR,grabRelR] = checkLevergrab(tapDur,postStimDur,levGrabR{iTrials},levReleaseR{iTrials},1/sRate);
    grabR(grabR < -preStimDur) = [];
    if ~isempty(grabR); grabRelR(grabRelR<grabR(1)) = []; end
    
    for iRegs = 0 : length(motorIdx)-1
        
        shiftGrabL = grabL - (mPreTime - (iRegs * 1/sRate));
        shiftGrabR = grabR - (mPreTime - (iRegs * 1/sRate));
        
        shiftGrabRelL = grabRelL - (mPreTime - (iRegs * 1/sRate));
        shiftGrabRelR = grabRelR - (mPreTime - (iRegs * 1/sRate));
        
        if iRegs == 0 %first regressor, find first grab / tap
            
            lGrabR{iTrials}(find(histcounts(shiftGrabL,-preStimDur:1/sRate:postStimDur),1),iRegs+1) = true;
            rGrabR{iTrials}(find(histcounts(shiftGrabR,-preStimDur:1/sRate:postStimDur),1),iRegs+1) = true;
            
        else
            
            regOut = leverEvents(iRegs, grabRelL, shiftGrabL, find(lGrabR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPreTime*sRate, 'pre'); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            lGrabR{iTrials}(regOut,iRegs+1) = true;
            
            regOut = leverEvents(iRegs, grabRelR, shiftGrabR, find(rGrabR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPreTime*sRate, 'pre'); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            rGrabR{iTrials}(regOut,iRegs+1) = true;
            
            regOut = leverEvents(iRegs, grabL, shiftGrabRelL, find(lGrabRelR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPreTime*sRate, ''); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            lGrabRelR{iTrials}(regOut,iRegs+1) = true;
            
            regOut = leverEvents(iRegs, grabR, shiftGrabRelR, find(rGrabRelR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPreTime*sRate, ''); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            rGrabRelR{iTrials}(regOut,iRegs+1) = true;           
            
        end
    end
        
    if gaussShift > 1
        % subsample regressors
        lGrabR{iTrials} = lGrabR{iTrials}(:,1:gaussShift:end);
        lGrabRelR{iTrials} = lGrabRelR{iTrials}(:,1:gaussShift:end);
        rGrabR{iTrials} = rGrabR{iTrials}(:,1:gaussShift:end);
        rGrabRelR{iTrials} = rGrabRelR{iTrials}(:,1:gaussShift:end);
    end
    
    %% pupil / whisk / nose / face regressors
    bhvFrameRate = round(1/mean(diff(pTime{bTrials(iTrials)}))); %framerate of face camera
    trialOn = bhv.TrialStartTime(iTrials) + (stimGrab(iTrials) - preStimDur);
    trialTime = pTime{bTrials(iTrials)} - trialOn;    
    idx = trialTime < trialDur; %don't use late frames
    trialTime = trialTime(idx);
    timeLeft = trialDur - trialTime(end); %check if there is missing time at the end of a trial
       
    if (timeLeft < trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
        addTime = trialTime(end) + (1/bhvFrameRate : 1/bhvFrameRate : timeLeft + 1/bhvFrameRate); %add some dummy times to make complete trial
        trialTime = [trialTime' addTime];
    end
    
    fastPupilR{iTrials} = Behavior_vidResamp(fPupil{bTrials(iTrials)}(idx), trialTime, sRate);
    fastPupilR{iTrials} = fastPupilR{iTrials}(end - frames + 1 : end);
        
    slowPupilR{iTrials} = Behavior_vidResamp(sPupil{bTrials(iTrials)}(idx), trialTime, sRate);
    slowPupilR{iTrials} = slowPupilR{iTrials}(end - frames + 1 : end);
        
    whiskR{iTrials} = Behavior_vidResamp(whisker{bTrials(iTrials)}(idx), trialTime, sRate);
    whiskR{iTrials} = smooth(whiskR{iTrials}(end - frames + 1 : end), 'rlowess');
        
    noseR{iTrials} = Behavior_vidResamp(nose{bTrials(iTrials)}(idx), trialTime, sRate);
    noseR{iTrials} = smooth(noseR{iTrials}(end - frames + 1 : end), 'rlowess');
    
    faceR{iTrials} = Behavior_vidResamp(faceM{bTrials(iTrials)}(idx), trialTime, sRate);
    faceR{iTrials} = smooth(faceR{iTrials}(end - frames + 1 : end), 'rlowess');
    
    %% body regressors
    bhvFrameRate = round(1/mean(diff(bTime{bTrials(iTrials)}))); %framerate of body camera
    trialTime = bTime{bTrials(iTrials)} - trialOn;
    idx = trialTime < trialDur; %don't use late frames
    trialTime = trialTime(idx);
    timeLeft = trialDur - trialTime(end); %check if there is missing time at the end of a trial
    
    if (timeLeft < trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
        addTime = trialTime(end) + (1/bhvFrameRate : 1/bhvFrameRate : timeLeft + 1/bhvFrameRate); %add some dummy times to make complete trial
        trialTime = [trialTime' addTime];
    end
    
    bodyR{iTrials} = Behavior_vidResamp(bodyM{bTrials(iTrials)}(idx), trialTime, sRate);
    bodyR{iTrials} = smooth(bodyR{iTrials}(end - frames + 1 : end), 'rlowess');        
        
    %% piezo sensor information
    if exist([cPath 'Analog_'  num2str(trials(iTrials)) '.dat'],'file') ~= 2  %check if files exists on hdd and pull from server otherwise
        cFile = dir([sPath 'Analog_'  num2str(trials(iTrials)) '.dat']);
        copyfile([sPath 'Analog_'  num2str(trials(iTrials)) '.dat'],[cPath 'Analog_'  num2str(trials(iTrials)) '.dat']);
    end
    
    [~,Analog] = Widefield_LoadData([cPath 'Analog_'  num2str(trials(iTrials)) '.dat'],'Analog'); %load analog data    
    stimOn = find(diff(double(Analog(stimLine,:)) > 1500) == 1); %find stimulus onset in current trial
    
    Analog(1,round(stimOn + ((postStimDur-stimTime(iTrials)) * 1000) - 1)) = 0; %make sure there are enough datapoints in analog signal
    temp = Analog(piezoLine,round(stimOn - ((preStimDur + stimTime(1)) * 1000)) : round(stimOn + ((postStimDur - stimTime(1))* 1000) - 1)); % data from piezo sensor. Should encode animals hindlimb motion.
    temp = smooth(double(temp), sRate*5, 'lowess')'; %do some smoothing
    temp = [repmat(temp(1),1,1000) temp repmat(temp(end),1,1000)]; %add some padding on both sides to avoid edge effects when resampling
    temp = resample(double(temp), sRate, 1000); %resample to imaging rate
    piezoR{iTrials} = temp(sRate + 1 : end - sRate)'; %remove padds again

    temp = abs(hilbert(diff(piezoR{iTrials})));
    piezoMoveR{iTrials} = [temp(1); temp]; %keep differential motion signal
    clear temp

    % give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,trialCnt);
        toc
    end
end

%% rebuild analog motor regressors to get proper design matrices
[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,whiskR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
temp = [single(traceOut) cat(1,dMat{:})]; %rebuild continuous format
[dMat, ~] = Widefield_analogToDesign(double(cat(1,whiskR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
whiskR = [temp cat(1,dMat{:})]; %add high amplitude movements separately

[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,noseR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
temp = [single(traceOut) cat(1,dMat{:})]; %rebuild continuous format
[dMat, ~] = Widefield_analogToDesign(double(cat(1,noseR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
noseR = [temp cat(1,dMat{:})]; %add high amplitude movements separately

[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,piezoR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
piezoR1 = [traceOut cat(1,dMat{:})]; %rebuild continuous format
[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,piezoMoveR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
temp = [traceOut cat(1,dMat{:})]; %rebuild continuous format
[dMat, ~] = Widefield_analogToDesign(double(cat(1,piezoMoveR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
piezoR = [piezoR1 temp cat(1,dMat{:})]; %add high amplitude movements separately

[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,faceR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
temp = [single(traceOut) cat(1,dMat{:})]; %rebuild continuous format
[dMat, ~] = Widefield_analogToDesign(double(cat(1,faceR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
faceR = [temp cat(1,dMat{:})]; %add high amplitude movements separately

[dMat, traceOut] = Widefield_analogToDesign(double(cat(1,bodyR{:})), 0.5, trialCnt, sRate, sRate, motorIdx, gaussShift);
temp = [single(traceOut) cat(1,dMat{:})]; %rebuild continuous format
[dMat, ~] = Widefield_analogToDesign(double(cat(1,bodyR{:})), 2, trialCnt, sRate, sRate, motorIdx, gaussShift);
bodyR = [temp cat(1,dMat{:})]; %add high amplitude movements separately
clear piezoR1 piezoR2 dMat traceOut temp

%% re-align behavioral video data and Vc to lever grab instead of stimulus onset
iiSpikeFrames = findInterictalSpikes(U, Vc); %find interictal spikes
Vc = interpOverInterictal(Vc, iiSpikeFrames); %interpolate over interictal spikes

V1 = reshape(V1,205,[],bhvDimCnt); %get to trial format
V2 = reshape(V2,205,[],bhvDimCnt); %get to trial format
vidR = V1(:,bTrials,:); clear V1 %get correct trials from behavioral video data.
moveR = V2(:,bTrials,:); clear V2 %get correct trials from behavioral video data.

% re-align video data
temp1 = NaN(dims,frames,trialCnt);
temp2 = NaN(frames,trialCnt,bhvDimCnt);
temp3 = NaN(frames,trialCnt,bhvDimCnt);
for x = 1 : size(vidR,2)
    try
        temp1(:,:,x) = Vc(:,(91 - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (91 - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x);
        temp2(:,x,:) = vidR((91 - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (91 - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x,:);
        temp3(:,x,:) = moveR((91 - ceil(stimTime(x) / (1/sRate))) - (preStimDur * sRate) : (91 - ceil(stimTime(x) / (1/sRate))) + (postStimDur * sRate) - 1,x,:);
    catch
        fprintf(1,'Could not align trial %d. Relative stim time: %fs\n', x, stimTime(x));
    end
end
Vc = reshape(temp1,dims,[]); clear temp1
vidR = reshape(temp2,[],bhvDimCnt); clear temp2
moveR = reshape(temp3,[],bhvDimCnt); clear temp2

%% reshape regressors, make design matrix and indices for regressors that are used for the model
timeR = repmat(logical(diag(ones(1,frames))),trialCnt,1); %time regressor
timeR = timeR(:,1:end-4);

lGrabR = cat(1,lGrabR{:});
lGrabRelR = cat(1,lGrabRelR{:});
rGrabR = cat(1,rGrabR{:});
rGrabRelR = cat(1,rGrabRelR{:});

lLickR = cat(1,lLickR{:});
rLickR = cat(1,rLickR{:});
leverInR = cat(1,leverInR{:});

lVisStimR = cat(1,lVisStimR{:});
rVisStimR = cat(1,rVisStimR{:});
lAudStimR = cat(1,lAudStimR{:});
rAudStimR = cat(1,rAudStimR{:});
spoutR = cat(1,spoutR{:});
spoutOutR = cat(1,spoutOutR{:});

visRewardR = cat(1,visRewardR{:});
audRewardR = cat(1,audRewardR{:});
prevRewardR = cat(1,prevRewardR{:});

ChoiceR = cat(1,ChoiceR{:});

prevChoiceR = cat(1,prevChoiceR{:});
prevModR = cat(1,prevModR{:});

waterR = cat(1,waterR{:});

fastPupilR = cat(1,fastPupilR{:});
fastPupilR(~isnan(fastPupilR(:,1)),:) = zscore(fastPupilR(~isnan(fastPupilR(:,1)),:));

slowPupilR = cat(1,slowPupilR{:});
slowPupilR(~isnan(slowPupilR(:,1)),:) = zscore(slowPupilR(~isnan(slowPupilR(:,1)),:));

if ~any(isnan(fastPupilR(:,1)) == isnan(vidR(:,1)))
    error('Pupil and bhv video do not agree with one another. Maybe something wrong with timing ?')
end

%% create full design matrix
fullR = [timeR ChoiceR visRewardR audRewardR lGrabR lGrabRelR rGrabR rGrabRelR lLickR rLickR lVisStimR rVisStimR lAudStimR rAudStimR ...
    prevRewardR prevChoiceR prevModR waterR piezoR whiskR noseR fastPupilR slowPupilR faceR bodyR moveR vidR];

% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
recLabels = {
    'time' 'Choice' 'visReward' 'audReward' 'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick' 'lVisStim' 'rVisStim' ...
    'lAudStim' 'rAudStim' 'prevReward' 'prevChoice' 'prevMod' 'water' 'piezo' 'whisk' 'nose' 'fastPupil' 'slowPupil' 'face' 'body' 'Move' 'bhvVideo'};

%index to reconstruct different response kernels
recIdx = [
    ones(1,size(timeR,2))*find(ismember(recLabels,'time')) ...
    ones(1,size(ChoiceR,2))*find(ismember(recLabels,'Choice')) ...
    ones(1,size(visRewardR,2))*find(ismember(recLabels,'visReward')) ...
    ones(1,size(audRewardR,2))*find(ismember(recLabels,'audReward')) ...
    ones(1,size(lGrabR,2))*find(ismember(recLabels,'lGrab')) ...
    ones(1,size(lGrabRelR,2))*find(ismember(recLabels,'lGrabRel')) ...
    ones(1,size(rGrabR,2))*find(ismember(recLabels,'rGrab')) ...
    ones(1,size(rGrabRelR,2))*find(ismember(recLabels,'rGrabRel')) ...
    ones(1,size(lLickR,2))*find(ismember(recLabels,'lLick')) ...
    ones(1,size(rLickR,2))*find(ismember(recLabels,'rLick')) ...
    ones(1,size(lVisStimR,2))*find(ismember(recLabels,'lVisStim')) ...
    ones(1,size(rVisStimR,2))*find(ismember(recLabels,'rVisStim')) ...
    ones(1,size(lAudStimR,2))*find(ismember(recLabels,'lAudStim')) ...
    ones(1,size(rAudStimR,2))*find(ismember(recLabels,'rAudStim')) ...
    ones(1,size(prevRewardR,2))*find(ismember(recLabels,'prevReward')) ...
    ones(1,size(prevChoiceR,2))*find(ismember(recLabels,'prevChoice')) ...
    ones(1,size(prevModR,2))*find(ismember(recLabels,'prevMod')) ...
    ones(1,size(waterR,2))*find(ismember(recLabels,'water')) ...
    ones(1,size(piezoR,2))*find(ismember(recLabels,'piezo')) ...
    ones(1,size(whiskR,2))*find(ismember(recLabels,'whisk')) ...
    ones(1,size(noseR,2))*find(ismember(recLabels,'nose')) ...
    ones(1,size(fastPupilR,2))*find(ismember(recLabels,'fastPupil')) ...
    ones(1,size(slowPupilR,2))*find(ismember(recLabels,'slowPupil')) ...
    ones(1,size(faceR,2))*find(ismember(recLabels,'face')) ...
    ones(1,size(bodyR,2))*find(ismember(recLabels,'body')) ...
    ones(1,size(moveR,2))*find(ismember(recLabels,'Move')) ...
    ones(1,size(vidR,2))*find(ismember(recLabels,'bhvVideo'))];

trialIdx = isnan(mean(fullR,2)); %don't use first trial or trials that failed to contain behavioral video data
fprintf(1, 'Rejected %d/%d trials for NaN entries in regressors\n', sum(trialIdx)/frames,trialCnt);
fullR(trialIdx,:) = []; %clear bad trials

idx = nansum(abs(fullR)) < 10; %reject regressors that are too sparse
fullR(:,idx) = []; %clear empty regressors
fprintf(1, 'Rejected %d/%d empty regressors\n', sum(idx),length(idx));

%% run QR and check for rank-defficiency
[~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize design matrix
figure; plot(abs(diag(fullQRR))); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
    error('Design matrix is rank-defficient')
end

%% save modified Vc
Vc(:,trialIdx) = []; %clear bad trials
save([cPath 'interpVc.mat'], 'Vc', 'iiSpikeFrames','frames');

%% apply gaussian filter to design matrix if using sub-sampling
if gaussShift > 1
    [a,b] = size(fullR);
    
    % find non-continous regressors (contain values different from -1, 0 or 1)
    temp = false(size(fullR));
    temp(fullR(:) ~= 0 & fullR(:) ~= 1 & fullR(:) ~= -1 & ~isnan(fullR(:))) = true;
    regIdx = nanmean(temp) == 0; %index for non-continous regressors
    
    % do gaussian convolution. perform trialwise to avoid overlap across trials.
    trialCnt = a/frames;
    fullR = reshape(fullR,frames,trialCnt,b);
    for iTrials = 1:trialCnt
        fullR(:,iTrials,regIdx) = smoothCol(squeeze(fullR(:,iTrials,regIdx)),gaussShift*2,'gauss');
    end
    fullR = reshape(fullR,a,b);
end

%% clear individual regressors
clear stimR lGrabR lGrabRelR rGrabR rGrabRelR waterR lLickR rLickR ...
    lVisStimR rVisStimR lAudStimR rAudStimR visRewardR audRewardR prevRewardR visChoiceR audChoiceR ...
    prevChoiceR prevModR fastPupilR moveR piezoR whiskR noseR faceR bodyR

%% run ridge regression in low-D
%run model. Zero-mean without intercept. no QR.
[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for original video, zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'orgdimBeta.mat'], 'dimBeta', 'ridgeVals')
save([cPath filesep 'orgregData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels','gaussShift','fullQRR','-v7.3');
Behavior_betaRebuild(cPath, 'org'); % rebuild video regressors by projecting beta weights for each wiedfield dimensions back on the behavioral video data

%% orthogonalize some regressors for clarity
% orthogonalize motor regressors against pupil to reduce state correlations.
pInd = ismember(recIdx(~idx), find(ismember(recLabels,{'fastPupil', 'slowPupil'})));
lInd = ismember(recIdx(~idx), find(ismember(recLabels,{'lLick', 'rLick'})));
wInd = ismember(recIdx(~idx), find(ismember(recLabels,'whisk')));
nInd = ismember(recIdx(~idx), find(ismember(recLabels,'nose')));
fInd = ismember(recIdx(~idx), find(ismember(recLabels,'face')));

smallR = [fullR(:,pInd) fullR(:,lInd) fullR(:,wInd) fullR(:,nInd) fullR(:,fInd)];
[Q, ~] = qr(smallR,0); clear smallR %orthogonalize motor from pupil
fullR(:,lInd) = Q(:,sum(pInd) + 1 : sum(pInd) + sum(lInd));
fullR(:,wInd) = Q(:,sum(pInd) + sum(lInd) + 1 : sum(pInd) + sum(lInd) + sum(wInd));
fullR(:,nInd) = Q(:,sum(pInd) + sum(lInd) + sum(wInd) + 1 : sum(pInd) + sum(lInd) + sum(wInd) + sum(nInd));
fullR(:,fInd) = Q(:,sum(pInd) + sum(lInd) + sum(wInd) + sum(nInd) + 1 : end);

% orthogonlize piezo/body from pupils
piInd = ismember(recIdx(~idx), find(ismember(recLabels,{'piezo'})));
smallR = [fullR(:,pInd) fullR(:,piInd)];
[Q, ~] = qr(smallR,0); clear smallR %orthogonalize piezo from pupil
fullR(:,piInd) = Q(:,sum(pInd) + 1 : sum(piInd) + sum(pInd));

% orthogonalize video against other motor regressors and spout movement
motorLabels = {'lGrab', 'rGrab', 'lGrabRel', 'rGrabRel', 'lLick', 'rLick', 'piezo', 'whisk', 'nose', 'fastPupil', 'slowPupil', 'face', 'body', 'Move', 'bhvVideo'};
cInd = ismember(recIdx(~idx), find(ismember(recLabels,motorLabels)));
smallR = [leverInR(~trialIdx,:) spoutR(~trialIdx,:) spoutOutR(~trialIdx,:) fullR(:,cInd)];

[Q, redQRR] = qr(smallR,0); clear smallR %orthogonalize video against motor regressors
vidIdx = recIdx(~idx) == find(ismember(recLabels,'bhvVideo'));
moveIdx = recIdx(~idx) == find(ismember(recLabels,'Move'));

% transfer orthogonolized video regressors back to design matrix
temp = size(leverInR,2) + size(spoutR,2) + size(spoutOutR,2); %length of non-motor regressors
fullR(:,vidIdx) = Q(:,[false(1,temp) vidIdx(cInd)]);
fullR(:,moveIdx) = Q(:,[false(1,temp) moveIdx(cInd)]);

%% run model with orthogonalized video. Zero-mean without intercept.
[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'dimBeta.mat'], 'dimBeta', 'ridgeVals')
save([cPath filesep 'regData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels','gaussShift','redQRR','-v7.3');
Behavior_betaRebuild(cPath); % rebuild video regressors by projecting beta weights for each wiedfield dimensions back on the behavioral video data

%% run video only model. Zero-mean without intercept.
fullR = vidR(~trialIdx,:);
idx = false(1,size(vidR,2));
recLabels = {'bhvVideo'};
recIdx = ones(1,size(vidR,2));
[ridgeVals, dimBeta] = ridgeMML(Vc', fullR, true); %get ridge penalties and beta weights.
fprintf('Mean ridge penalty for video-only zero-mean model: %f\n', mean(ridgeVals));
save([cPath 'vidOnlydimBeta.mat'], 'dimBeta', 'ridgeVals')
save([cPath filesep 'vidOnlyregData.mat'], 'fullR', 'idx', 'trialIdx', 'recIdx', 'recLabels','gaussShift','-v7.3');
Behavior_betaRebuild(cPath, 'vidOnly'); % rebuild video regressors by projecting beta weights for each wiedfield dimensions back on the behavioral video data


%% nested functions
function regOut = leverEvents(cReg, times, shiftTimes, lastTimes, timeFrame, preTimeShift, dMode) %code to compute events for current lever regressor

regOut = [];
if cReg == 0 %first regressor, find first grab
    regOut = find(histcounts(shiftTimes,timeFrame),1);
    
elseif cReg < preTimeShift && strcmpi(dMode,'pre') %negative shift regressors - only for grab or tap
    regOut = lastTimes + 1; %find last event and shift one forward
    if isempty(regOut)
        regOut = find(histcounts(shiftTimes,timeFrame),1); %check if new event can be found if none is present so far
    end
    
elseif cReg == preTimeShift %this is the zero-lag regressor. use all available events
    regOut = find(histcounts(shiftTimes,timeFrame)); %times and shifttimes should be the same here.
    
elseif cReg > preTimeShift %for positive shift, use the last regressor but eliminate event if there is overlap with zero-lag regressor.
    regOut = lastTimes + 1; %find last event and shift one forward
    if ~isempty(regOut)
        regOut(ismember(regOut, find(histcounts(times,timeFrame)))) = []; %remove events that overlap with the given zero-lag regressor. Can also be from a different event type like handle release.
    end
end



function [grabOn,grabRel,tapOn,tapRel] = checkLevergrab(tapDur,postStimDur,grabs,release,minTime)

grabOn = [];
grabRel = [];
tapOn = [];
tapRel = [];

if ~isempty(grabs)
    Cnt = 0;
    grabCnt = 0;
    tapCnt = 0;
    
    while true
        Cnt = Cnt +1;
        if length(grabs) < Cnt %no more grabs
            break;
        end
        
        idx = find(release > grabs(Cnt),1);
        if ~isempty(idx) %if there is a release that follows current grab
            cGrab = release(idx) - grabs(Cnt); %duration of current grab
        else
            cGrab = postStimDur - grabs(Cnt);
        end
        
        %% more grabs available - check start time of next grab and merge if they are very close in time
        if length(grabs) > Cnt && ~isempty(idx)
            while (grabs(Cnt+1) - release(idx)) <= minTime %time diff between grabs is less than minimum duration. Merge with next grab.
                
                release(idx) = []; %delete current release
                grabs(Cnt+1) = []; %delete next grab
                
                idx = find(release > grabs(Cnt),1); %next release
                if ~isempty(idx) %if there is a release that follows current grab
                    cGrab = release(idx) - grabs(Cnt); %duration of current grab
                else
                    cGrab = postStimDur - grabs(Cnt);
                end
                
                if length(grabs) <= Cnt %no more grabs
                    break
                end
                
            end
        end
        
        %% check if current grab is grab or tap
        if cGrab <= tapDur
            tapCnt = tapCnt + 1;
            tapOn(tapCnt) = grabs(Cnt);
            if ~isempty(idx)
                tapRel(tapCnt) = idx;
            end
        else
            grabCnt = grabCnt + 1;
            grabOn(grabCnt) = grabs(Cnt);
            if ~isempty(idx)
                grabRel(grabCnt) = release(idx);
            end
        end
        
        if isempty(idx) || length(grabs) <= Cnt %no more grabs/releases
            break;
        end
    end
end