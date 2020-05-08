function BpodImager_squareRegressModel(cPath,Animal,Rec)

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

%general variables
Paradigm = 'SpatialDisc';
cPath = [cPath Animal filesep Paradigm filesep Rec filesep]; %Widefield data path
sPath = ['U:\space_managed_data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
sRate = 40;         % Sampling rate of imaging in Hz

% trial structure variables
preStimDur = 2;                 % Duration of trial before stimulus onset in seconds
postStimDur = 2;                % Duration of trial after stimulus onset in seconds
trialSegments = ones(1,4);  % Duration of  different segments in each trial in seconds. This refers to baseline, lever grab, stimulus and decision phase and usually there are 1s each.

% ridge regression variables
ridgeCycles = 5;    %number of cycles for ridge estimation
ridgeSteps = 10;    %number of steps for first cycle
ridgeFolds = 2;     %number of folds for cross validation
maxRidge = 50;       %range of values around ridgeVal that should be tested
ridgeVal = 200;       %initial estimate for ridge parameter

%other variables
mPretime = 1;       % precede self-initiated events in case a neuron predicts animal motor action
tapDur = 0.25;      % minimum time of lever contact, required to count as a proper grab. Will be counted as 'tap' otherwise.
bhvDimCnt = 100;    % number of dimensions from behavioral videos that are used as regressors.
smth = 5;           % filter length when smoothing video data
smth = smth-1+mod(smth,2); % ensure kernel length is odd

newRun = true; %flag to indicate that all values should be recomputed (including ridge regression and correlation maps)

%% load behavior data
if exist([cPath 'BehaviorVideo' filesep 'SVD_Cam 1.mat'],'file') ~= 2 %check if svd behavior exists on hdd and pull from server otherwise
    if exist([cPath 'BehaviorVideo' filesep]) ~= 2
        mkdir([cPath 'BehaviorVideo' filesep]);
    end
    copyfile([sPath 'BehaviorVideo' filesep 'SVD_Cam 1.mat'],[cPath 'BehaviorVideo' filesep 'SVD_Cam 1.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'SVD_Cam 2.mat'],[cPath 'BehaviorVideo' filesep 'SVD_Cam 2.mat']);
end

load([cPath 'BehaviorVideo' filesep 'SVD_Cam 1.mat']);
% V1 = zscore(V')./10; 
V1 = V'; 
U1 = U;
frameTimes1 = totalFrameTimes;

load([cPath 'BehaviorVideo' filesep 'SVD_Cam 2.mat']);
% V2 = zscore(V')./10;
V2 = V';
U2 = U;
frameTimes2 = totalFrameTimes;

%% load data
load([cPath 'mask.mat'])
load([cPath 'Vc.mat'])

bhvFile = dir([cPath filesep Animal '_' Paradigm '*.mat']);
load([cPath bhvFile(1).name]); %load behavior data

% ensure there are not too many trials in Vc
ind = trials > SessionData.nTrials;
trials(ind) = [];
Vc(:,:,ind) = [];
[dims, frames, trialCnt] = size(Vc);

bhv = selectBehaviorTrials(SessionData,trials); %only use completed trials that are in the Vc dataset

% %% get time point of stimulus onset from analog data
% stimOn = BpodImager_recoverStimOn(cPath,sPath,trials);

%% find events in BPod time.
% All timestamps are relative to stimulus onset event to synchronize to imaging data later

% pre-allocate vectors
lickL = cell(1,trialCnt);
lickR = cell(1,trialCnt);
leverIn = NaN(1,trialCnt);
visStim = NaN(1,trialCnt);
audStim = NaN(1,trialCnt);
levGrabL = cell(1,trialCnt);
levGrabR = cell(1,trialCnt);
levReleaseL = cell(1,trialCnt);
levReleaseR = cell(1,trialCnt);
water = NaN(1,trialCnt);

for iTrials = 1:trialCnt
    
    if bhv.StimType(iTrials) == 1 || bhv.StimType(iTrials) == 3 %visual or mixed stimulus
        try
            visStim(iTrials) = bhv.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
        end
    end
    if bhv.StimType(iTrials) == 2 || bhv.StimType(iTrials) == 3 %auditory or mixed stimulus
        try
            audStim(iTrials) = bhv.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
        end
    end
    stimTime = nanmean([visStim(iTrials) audStim(iTrials)]);
    
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port1In') %check for licks
        lickL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port1In;
        lickL{iTrials}(lickL{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickL{iTrials} = lickL{iTrials} - stimTime;
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port3In') %check for right licks
        lickR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port3In;
        lickR{iTrials}(lickR{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickR{iTrials} = lickR{iTrials} - stimTime;
    end
    
    leverIn(iTrials) = min(bhv.RawEvents.Trial{iTrials}.States.Reset(:)) - stimTime; %first reset state causes lever to move in
    
    leverTimes = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
    
    stimGrab(iTrials) = leverTimes(find(leverTimes == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1) - stimTime; %find start of lever state that triggered stimulus onset
    
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2High') %check for left grabs
        levGrabL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2High - stimTime;
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1High') %check for right grabs
        levGrabR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1High - stimTime;
    end
    
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1Low') %check for left release
        levReleaseL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1Low - stimTime;
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2Low') %check for right release
        levReleaseR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2Low - stimTime;
    end
    
    if ~isnan(bhv.RawEvents.Trial{iTrials}.States.Reward(1)) %check for reward state
        water(iTrials) = bhv.RawEvents.Trial{iTrials}.States.Reward(1) - stimTime;
    end
end

bhv.ResponseSide(bhv.ResponseSide == 1) = -1; %normalize response side between -1 and 1
bhv.ResponseSide(bhv.ResponseSide == 2) = 1; %normalize response side between -1 and 1

%% build regressors - create design matrix based on event times
% spoutMove = round(bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1) - bhv.RawEvents.Trial{iTrials}.Events.Wire3High(1),3); %round to ms
lickWindow = (postStimDur - (1 - mPretime)) * sRate; %time window to be covered by lick regressors

%create basic time regressors
timeR = logical(diag(ones(1,frames)));

% allocate cells for other regressors
leverAlignR = cell(1,trialCnt);
stimAlignR = cell(1,trialCnt);
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

visRewardR = cell(1,trialCnt);
visPrevRewardR = cell(1,trialCnt);
audRewardR = cell(1,trialCnt);
audPrevRewardR = cell(1,trialCnt);
mixRewardR = cell(1,trialCnt);
mixPrevRewardR = cell(1,trialCnt);

choiceR = cell(1,trialCnt);
prevChoiceR = cell(1,trialCnt);

waterR = cell(1,trialCnt);

pupilR = cell(1,trialCnt);
faceR = cell(1,trialCnt);
bodyR = cell(1,trialCnt);

tic
for iTrials = 1:trialCnt
    % baseline/lever (leverAlignR) and stim/decision (stimAlignR) regressors - these are required to enable variable stimulus onset times
    
    % time regressor, aligned to start of the (last) lever period. Zero lag regressor at frames / 4.
    offset = (frames / 4) - find(histcounts(stimGrab(iTrials),-preStimDur:1/sRate:postStimDur));
    leverAlignR{iTrials} = false(frames,frames);
    leverAlignR{iTrials}(:,offset+1:end) = timeR(:,1:(frames) - offset);
    leverAlignR{iTrials}((frames/2)+1:end,:) = false; %should not extent into stimulus period at frames / 2
    
    stimAlignR{iTrials} = timeR(:,(size(timeR,2)/2)+1:end); %time regressor, aligned to and starting at stim onset.
    
    %% vis/aud stim
    %regressors cover the last 2s after stimOn
    if ~isnan(visStim(iTrials)) && ~isempty(visStim(iTrials))
        if bhv.CorrectSide(iTrials) == 1
            lVisStimR{iTrials} = timeR(:,(size(timeR,2)/2)+1: (size(timeR,2)/2) + (postStimDur * sRate));
            rVisStimR{iTrials} = false(size(timeR,1), postStimDur * sRate);
        else
            lVisStimR{iTrials} = false(size(timeR,1), postStimDur * sRate);
            rVisStimR{iTrials} = timeR(:,(size(timeR,2)/2)+1: (size(timeR,2)/2) + (postStimDur * sRate));
        end
    else
        lVisStimR{iTrials} = false(size(timeR,1), postStimDur * sRate);
        rVisStimR{iTrials} = false(size(timeR,1), postStimDur * sRate);
    end
    
    if ~isnan(audStim(iTrials)) && ~isempty(audStim(iTrials))
        if bhv.CorrectSide(iTrials) == 1
            lAudStimR{iTrials} = timeR(:,(size(timeR,2)/2)+1: (size(timeR,2)/2) + (postStimDur * sRate));
            rAudStimR{iTrials} = false(size(timeR,1), postStimDur * sRate);
        else
            lAudStimR{iTrials} = false(size(timeR,1), postStimDur * sRate);
            rAudStimR{iTrials} = timeR(:,(size(timeR,2)/2)+1: (size(timeR,2)/2) + (postStimDur * sRate));
        end
    else
        lAudStimR{iTrials} = false(size(timeR,1), postStimDur * sRate);
        rAudStimR{iTrials} = false(size(timeR,1), postStimDur * sRate);
    end
    
    %% lick regressors
    lLickR{iTrials} = false(frames, lickWindow + 1);
    rLickR{iTrials} = false(frames, lickWindow + 1);
    
    for iRegs = 0 : lickWindow
        licks = lickL{iTrials} - (mPretime - (iRegs * 1/sRate));
        lLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
        
        licks = lickR{iTrials} - (mPretime - (iRegs * 1/sRate));
        rLickR{iTrials}(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
    end
    
    %% lever in
    leverInR{iTrials} = false(size(timeR,1),size(timeR,2));
    leverShift = (size(timeR,2)/2) + round(leverIn(iTrials) * sRate); %timepoint in frames were lever moved in, relative to stimOnset
    
    if ~isnan(leverShift)
        if leverShift > 0 %lever moved in during the recorded trial
            leverInR{iTrials}(:, 1: size(timeR,2) - leverShift) = timeR(:, leverShift+1:end);
        else %lever was present before data was recorded
            leverInR{iTrials}(:, abs(leverShift) + 1 : end) = timeR(:,  1: size(timeR,2) + leverShift);
        end
    end

    %% choice and reward
    visRewardR{iTrials} = (timeR .* (bhv.StimType(iTrials)  == 1)) * (-1 + (2 * bhv.Rewarded(iTrials))); %visual trial, 1 for reward, -1 for no reward
    audRewardR{iTrials} = (timeR .* (bhv.StimType(iTrials)  == 2)) * (-1 + (2 * bhv.Rewarded(iTrials))); %audio trial, 1 for reward, -1 for no reward
    mixRewardR{iTrials} = (timeR .* (bhv.StimType(iTrials)  == 3)) * (-1 + (2 * bhv.Rewarded(iTrials))); %audio trial, 1 for reward, -1 for no reward
    
    choiceR{iTrials} = single(timeR .* bhv.ResponseSide(iTrials));

    if iTrials == 1 %don't use first trial
        visPrevRewardR{iTrials} = NaN(size(timeR));
        audPrevRewardR{iTrials} = NaN(size(timeR));
        mixPrevRewardR{iTrials} = NaN(size(timeR));
        prevChoiceR{iTrials} = NaN(size(timeR));
    else %for all subsequent trials, use results of last trial to compute previous reward/choice
        visPrevRewardR{iTrials} = visRewardR{iTrials-1};
        audPrevRewardR{iTrials} = audRewardR{iTrials-1};
        mixPrevRewardR{iTrials} = mixRewardR{iTrials-1};
        prevChoiceR{iTrials} = choiceR{iTrials-1};        
    end
    
    %determine timepoint of reward given
    waterR{iTrials} = false(size(timeR,1), sRate);

    if ~isnan(water(iTrials)) && ~isempty(water(iTrials))
        waterOn = (size(timeR,2)/2) + round(water(iTrials) * sRate); %timepoint in frames when reward was given
        waterR{iTrials}(:, 1: size(timeR,2) - waterOn) = timeR(:, waterOn+1:end);
    end
    
    %% lever grab / release
    lGrabR{iTrials} = false(frames, frames + 1);
    lGrabRelR{iTrials} = false(frames, frames + 1);
    
    rGrabR{iTrials} = false(frames, frames + 1);
    rGrabRelR{iTrials} = false(frames, frames + 1);
    
    [grabL,grabRelL] = checkLevergrab(tapDur,postStimDur,levGrabL{iTrials},levReleaseL{iTrials});
    grabL(grabL < -preStimDur) = []; %ensure there are no grab events too early in the trial
    if ~isempty(grabL);grabRelL(grabRelL<grabL(1)) = []; end %make sure there are no release events before first grab
    
    [grabR,grabRelR] = checkLevergrab(tapDur,postStimDur,levGrabR{iTrials},levReleaseR{iTrials});
    grabR(grabR < -preStimDur) = [];
    if ~isempty(grabR); grabRelR(grabRelR<grabR(1)) = []; end
    
    for iRegs = 0:size(timeR,2)
        
        shiftGrabL = grabL - (mPretime - (iRegs * 1/sRate));
        shiftGrabR = grabR - (mPretime - (iRegs * 1/sRate));
        
        shiftGrabRelL = grabRelL - (mPretime - (iRegs * 1/sRate));
        shiftGrabRelR = grabRelR - (mPretime - (iRegs * 1/sRate));
        
        if iRegs == 0 %first regressor, find first grab / tap
            
            lGrabR{iTrials}(find(histcounts(shiftGrabL,-preStimDur:1/sRate:postStimDur),1),iRegs+1) = true;
            rGrabR{iTrials}(find(histcounts(shiftGrabR,-preStimDur:1/sRate:postStimDur),1),iRegs+1) = true;
            
        else
            
            regOut = leverEvents(iRegs, grabRelL, shiftGrabL, find(lGrabR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPretime*sRate, 'pre'); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            lGrabR{iTrials}(regOut,iRegs+1) = true;
            
            regOut = leverEvents(iRegs, grabRelR, shiftGrabR, find(rGrabR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPretime*sRate, 'pre'); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            rGrabR{iTrials}(regOut,iRegs+1) = true;
            
            regOut = leverEvents(iRegs, grabL, shiftGrabRelL, find(lGrabRelR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPretime*sRate, ''); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            lGrabRelR{iTrials}(regOut,iRegs+1) = true;
            
            regOut = leverEvents(iRegs, grabR, shiftGrabRelR, find(rGrabRelR{iTrials}(:,iRegs)), -preStimDur:1/sRate:postStimDur, mPretime*sRate, ''); %code to compute events for current lever regressor
            regOut(regOut > frames) = [];
            rGrabRelR{iTrials}(regOut,iRegs+1) = true;           
            
        end
    end
    
    %% pupil/video regressors
    frameFile = dir([cPath 'BehaviorVideo' filesep '*frameTimes_'  num2str(trials(iTrials),'%04i') '*_1.mat']);
    if exist([cPath 'BehaviorVideo' filesep 'faceVars_' int2str(trials(iTrials)) '.mat'],'file') ~= 2 || isempty(frameFile)  %check if files exists on hdd and pull from server otherwise
        if exist([cPath 'BehaviorVideo' filesep]) ~= 2
            mkdir([cPath 'BehaviorVideo' filesep]);
        end
        frameFile = dir([sPath 'BehaviorVideo' filesep '*frameTimes_'  num2str(trials(iTrials),'%04i') '*_1.mat']);
        
        %get results from pupil analysis
        copyfile([sPath filesep 'BehaviorVideo' filesep 'faceVars_' int2str(trials(iTrials)) '.mat'],[cPath filesep 'BehaviorVideo' filesep 'faceVars_' int2str(trials(iTrials)) '.mat']);
        % copy frametimes from server to local
        copyfile([sPath filesep 'BehaviorVideo' filesep frameFile.name],[cPath filesep 'BehaviorVideo' filesep frameFile.name]);
    end
    
    load([cPath filesep 'BehaviorVideo' filesep 'faceVars_' int2str(trials(iTrials)) '.mat'])
    load([cPath filesep 'BehaviorVideo' filesep frameFile.name])
    
    %absolute trial onset time for bpod - use to synchronize video to behavioral events
    trialOn = bhv.TrialStartTime(iTrials) * 86400 + (nanmean([visStim(iTrials) audStim(iTrials)]) - preStimDur);
    frameTimes = frameTimes  * 86400 - trialOn;
    trialOn = find(frameTimes > 0,1); %trial on in frames
    bhvFrameRate = round(1/median(diff(frameTimes)));
    trialOn = round(trialOn * (sRate / bhvFrameRate)); %trial on in resampled frames
    
    %resample to match imaging data framerate
    try
        pupil = mean(eyeVars.axes,2); %pupil diameter
        idx = zscore(pupil) < -2; %index for blinks
        pupil(idx) = mean(pupil); %replace blinks with pupil average
        pupil = [repmat(pupil(1),21,1); pupil; repmat(pupil(end),21,1)]; %add some padding on both sides to avoid edge effects when resampling
        pupil = resample(pupil, sRate, bhvFrameRate); %resample to match imaging data framerate
        offset = ceil(21 * sRate / bhvFrameRate); %required offset to remove padding after resampling
        pupil = pupil(offset : end - offset); %remove padds
        pupil = smooth(pupil,'rlowess'); %do some smoothing
        pupilR{iTrials} = single(pupil(trialOn:trialOn + (size(timeR,1) - 1))); %only use trial-relevant frames
    catch
        pupilR{iTrials} = NaN(size(timeR,1),1,'single'); %can't use this trial
    end
    
    %get svd data from cam1 - this is usually the one that sees the animals face
    trialOn = bhv.TrialStartTime(iTrials) * 86400 + (nanmean([visStim(iTrials) audStim(iTrials)]) - preStimDur);
    cTimes = frameTimes1 * 86400 - trialOn;
    trialOn = find(cTimes > 0,1); %trial on in frames
    bhvFrameRate = round(1/median(diff(cTimes)));
    faceMotion = V1(trialOn : trialOn + ceil(size(timeR,1)*(bhvFrameRate / sRate)), 1:bhvDimCnt);
    try
        faceMotion = [repmat(faceMotion(1,:),21,1); faceMotion; repmat(faceMotion(end,:),21,1)]; %add some padding on both sides to avoid edge effects when resampling
        faceMotion = resample(double(faceMotion), sRate, bhvFrameRate); %resample to match imaging data framerate
        faceMotion = conv2(faceMotion',ones(1,smth)/smth,'same')'; %smooth trace with moving average of 'smth' points
        faceMotion = faceMotion(offset : end - offset,:); %remove padds
        faceR{iTrials} = single(faceMotion(1:size(timeR,1),:)); %only use trial-relevant frames
    catch
        faceR{iTrials} = NaN(size(timeR,1),bhvDimCnt, 'single'); %can't use this trial
    end
    
    %get svd data from cam2 - this is usually the one that sees the animal from the bottom
    trialOn = bhv.TrialStartTime(iTrials) * 86400 + (nanmean([visStim(iTrials) audStim(iTrials)]) - preStimDur);
    cTimes = frameTimes2 * 86400 - trialOn;
    trialOn = find(cTimes > 0,1); %trial on in frames
    bhvFrameRate = round(1/median(diff(cTimes)));
    bodyMotion = V2(trialOn : trialOn + ceil(size(timeR,1)*(bhvFrameRate / sRate)), 1:bhvDimCnt);
    
    try
        bodyMotion = [repmat(bodyMotion(1,:),21,1); bodyMotion; repmat(bodyMotion(end,:),21,1)]; %add some padding on both sides to avoid edge effects when resampling
        bodyMotion = resample(double(bodyMotion), sRate, bhvFrameRate); %resample to match imaging data framerate
        bodyMotion = conv2(bodyMotion',ones(1,smth)/smth,'same')'; %smooth trace with moving average of 'smth' points
        bodyMotion = bodyMotion(offset : end - offset,:); %remove padds
        bodyR{iTrials} = single(bodyMotion(1:size(timeR,1),:)); %only use trial-relevant frames
    catch
        bodyR{iTrials} = NaN(size(timeR,1), bhvDimCnt, 'single'); %can't use this trial
    end
    
    % give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,trialCnt);
        toc
    end
end

%% reshape regressors, make design matrix and indices for regressors that are used for the model
leverAlignR = cat(1,leverAlignR{:});
stimAlignR = cat(1,stimAlignR{:});

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

visRewardR = cat(1,visRewardR{:});
audRewardR = cat(1,audRewardR{:});
mixRewardR = cat(1,mixRewardR{:});
visPrevRewardR = cat(1,visPrevRewardR{:});
audPrevRewardR = cat(1,audPrevRewardR{:});
mixPrevRewardR = cat(1,mixPrevRewardR{:});

choiceR = cat(1,choiceR{:});
prevChoiceR = cat(1,prevChoiceR{:});

waterR = cat(1,waterR{:});

pupilR = cat(1,pupilR{:});
faceR = cat(1,faceR{:});
bodyR = cat(1,bodyR{:});

pupilR(~isnan(pupilR)) = zscore(pupilR(~isnan(pupilR)));
faceR(~isnan(faceR)) = zscore(faceR(~isnan(faceR)));
bodyR(~isnan(bodyR)) = zscore(bodyR(~isnan(bodyR)));

% full design matrix
fullR = [leverAlignR stimAlignR lGrabR lGrabRelR rGrabR rGrabRelR lLickR rLickR ...
    leverInR lVisStimR rVisStimR lAudStimR rAudStimR ...
    visRewardR audRewardR mixRewardR visPrevRewardR ... 
    audPrevRewardR mixPrevRewardR choiceR prevChoiceR waterR pupilR faceR bodyR];

idx = sum(abs(fullR)) < 10; %reject regressors that are too sparse
fullR(:,idx) = []; %clear empty regressors
fprintf(1, 'Rejected %d/%d empty regressors\n', sum(idx),length(idx));

%index to reconstruct different response kernels
recIdx = [ones(1,size(leverAlignR,2)) ones(1,size(stimAlignR,2))*2 ones(1,size(lGrabR,2))*3 ones(1,size(lGrabRelR,2))*4 ones(1,size(rGrabR,2))*5 ...
    ones(1,size(rGrabRelR,2))*6 ones(1,size(lLickR,2))*7 ones(1,size(rLickR,2))*8 ones(1,size(leverInR,2))*9 ...
    ones(1,size(lVisStimR,2))*10 ones(1,size(rVisStimR,2))*11 ones(1,size(lAudStimR,2))*12 ones(1,size(rAudStimR,2))*13 ...
    ones(1,size(visRewardR,2))*14 ones(1,size(audRewardR,2))*15 ones(1,size(mixRewardR,2))*16 ones(1,size(mixPrevRewardR,2))*17 ...
    ones(1,size(visPrevRewardR,2))*18 ones(1,size(audPrevRewardR,2))*19 ones(1,size(choiceR,2))*20 ones(1,size(prevChoiceR,2))*21 ...
    ones(1,size(waterR,2))*22 ones(1,size(pupilR,2))*23 ones(1,size(faceR,2))*24 ones(1,size(bodyR,2))*25]; 

recLabels = {
    'leverAlign' 'stimAlign' 'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick' 'leverIn' 'lVisStim' 'rVistStim' ...
    'lAudStim' 'rAudStim' 'visReward' 'audReward' 'mixReward' 'visPrevReward' 'audPrevReward' 'mixPrevReward' 'choice' ...
    'prevChoice' 'water' 'pupil' 'face' 'body'};

%clear individual regressors
clear leverAlignR stimAlignR lGrabR lGrabRelR rGrabR rGrabRelR waterR lLickR rLickR leverInR ...
      lVisStimR rVisStimR lAudStimR rAudStimR visRewardR audRewardR mixRewardR ... 
      visPrevRewardR audPrevRewardR mixPrevRewardR choiceR prevChoiceR pupilR faceR bodyR

trialIdx = isnan(mean(fullR,2)); %don't use first trial or trials that failed to contain behavioral video data
fullR(trialIdx,:) = []; %clear bad trials
fprintf(1, 'Rejected %d/%d trials for NaN regressors\n', sum(trialIdx) / frames,trialCnt);

% save some results
save([cPath filesep 'regData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels');

%% run ridge regression in low-D
U = arrayShrink(U,mask);
Vc = reshape(Vc,dims,[]);
Vc(:,trialIdx) = [];
Vc = bsxfun(@minus, Vc, mean(Vc, 2)); %make sure Vc is zero-mean

if exist([cPath 'ridgeTest.mat']) == 2 && ~newRun
    load([cPath 'ridgeTest.mat'])
else
    ridgeRange = (-maxRidge : maxRidge / (ridgeSteps-1) * 2 : maxRidge) /2 + ridgeVal; %range of values to be tested
    ridgeRange(ridgeRange <= 0) = [];
    
    if ~isunix; figure; hold on; end
    for iCycles = 1:ridgeCycles
        ridgeRMSE = NaN(ridgeFolds,length(ridgeRange));
        Cnt = 0;
        for iSteps = ridgeRange
            Cnt = Cnt + 1;
            foldCnt = floor(size(Vc,2) / ridgeFolds);
            
            for iFolds = 1:ridgeFolds
                
                cIdx = true(1,size(Vc,2));
                cIdx(((iFolds - 1)*foldCnt) + (1:foldCnt)) = false;
                [~, ridgeRMSE(iFolds,Cnt)] = ridgeFast(Vc', fullR, iSteps, cIdx); %perform ridge regression to get model error
                fprintf('iFold: %d/%d, iSteps: %d/%d, iCycles: %d/%d\n', iFolds, ridgeFolds, Cnt, length(ridgeRange), iCycles, ridgeCycles);
                
                if iFolds == ridgeFolds && iSteps == ridgeRange(end)
                    
                    ridgeTest{iCycles}{1} = ridgeRange;
                    ridgeTest{iCycles}{2} = ridgeRMSE;
                    if ~isunix; plot(ridgeRange, nanmean(ridgeRMSE),'o-','linewidth',2); axis square; end
                    
                    if iCycles ~= ridgeCycles
                        
                        [~,temp] = min(nanmean(ridgeRMSE)); %find ridge value that gave lowest error
                        if temp ~= size(ridgeRMSE,2)
                            maxRidge = ceil(maxRidge / 2); %decrease the range for next run
                        end
                        ridgeRange = (-maxRidge : maxRidge / (ridgeSteps-1) * 2 : maxRidge) /2 + ridgeRange(temp); %range of values to be tested
                        ridgeRange(ridgeRange <= 0) = [];
                        
                    end
                end
            end
        end
    end
    
    [~,temp] = min(nanmean(ridgeRMSE)); %find ridge value that gave lowest error
    ridgeVal = ridgeRange(temp); %use this value for regression on full dataset
    save([cPath 'ridgeTest.mat'], 'ridgeTest', 'ridgeVal')
end

if ~isunix
    plot([ridgeVal ridgeVal],get(gca,'YLim'),'k')
    text(ridgeVal+mean(get(gca,'XLim'))/50,mean(get(gca,'YLim')),'Ridge penalty')
    ylabel('Prediction RMSE'); xlabel('Ridge penalty');
    title(['Ridge penalty search; ' num2str(length(ridgeTest)) ' cycles - ' num2str(size(ridgeTest{1}{2},1)) 'xCrossVal, '])
end

%% run regression in low-D
[dimBeta, dimBetaCI] = ridgeFast(Vc', fullR, ridgeVal); %perform ridge regression to get beta weights and confidence intervals
beta = dimBeta * U'; %compute beta weight maps
betaCI = dimBetaCI * U'; %compute beta confidence maps
save([cPath 'dimBeta.mat'], 'dimBeta')
save([cPath 'dimBetaCI.mat'], 'dimBetaCI')

parBeta = cell(1,size(recLabels,2));
parBetaCI = cell(1,size(recLabels,2));
for iModels = 1:size(recLabels,2)
    cIdx = recIdx(~idx) == iModels;
    [parBeta{iModels}, parBetaCI{iModels}] = ridgeFast(Vc', fullR, ridgeVal, cIdx); %perform ridge regression to get model error
end
save([cPath 'parBeta.mat'], 'parBeta')
save([cPath 'parBetaCI.mat'], 'parBetaCI')

%% compute correlations in low-D
trialSegments = trialSegments * sRate;
if exist([cPath 'dimCorrFull.mat']) == 2 && ~newRun
    if ~isunix
        load([cPath 'dimCorrFull.mat'])
        figure;
        imagesc(arrayShrink(dimCorr.^2, mask, 'split')); axis square; colorbar; colormap jet
        title('R^2 - Dimensional full model')
    end
else
    Vm = (fullR * dimBeta)';
    dimCorrFull = matCorr(Vc, Vm, U); %compute correlation matrix between data and model estimate
    save([cPath 'dimCorrFull.mat'], 'dimCorrFull')
    
    for iModels = 1:length(parBeta)  %compute correlations for partial models
        parIdx = recIdx(~idx) == iModels;
        parCorrFull{iModels} = matCorr(Vc, (fullR(:,parIdx)*parBeta{iModels})', U); %compute correlation matrix between data and model estimate
    end
    save([cPath 'parCorrFull.mat'], 'parCorrFull', 'recLabels')
    
    dimCorrSeg = zeros(size(U,1),length(trialSegments));
    for iSegments = 1:length(trialSegments)
        % Get index for all frames that match the selected trial segments.
        % 1:40 baseline, 41:80 lever grab, 81:120 stimulus, 121:160 decision from each trial.
        a = ((0 : size(fullR,1) / sum(trialSegments) - 1) * sum(trialSegments));
        b = [0 cumsum(trialSegments)];
        cIdx = bsxfun(@plus,repmat(b(iSegments)+1:b(iSegments+1),length(a),1), a');
        cIdx = reshape(cIdx',1,[]);
        
        dimCorrSeg(:,iSegments) = matCorr(Vc(:,cIdx), Vm(:,cIdx), U); %compute correlation matrix between data and model estimate
        
        for iModels = 1:length(parBeta)  %compute correlations for partial models
            parIdx = (recIdx(~idx) == iModels); %index for current regressors
            parCorrSeg{iModels}(:,iSegments) = matCorr(Vc(:,cIdx), (fullR(cIdx,parIdx)*parBeta{iModels})', U); %compute correlation matrix between data and model estimate
        end
    end
    save([cPath 'dimCorrSeg.mat'], 'dimCorrSeg')
    save([cPath 'parCorrSeg.mat'], 'parCorrSeg', 'recLabels')
end

%% collect beta maps and reconstruct kernels
beta = arrayShrink(gather(beta)',mask,'split'); %spatially reconstruct beta maps
fullBeta = NaN(size(beta,1),size(beta,2),length(idx),'single');
fullBeta(:,:,~idx) = beta; %full matrix to include regressors that were rejected before

betaCI = arrayShrink(gather(betaCI)',mask,'split'); %spatially reconstruct beta confidence maps
fullBetaCI = NaN(size(betaCI,1),size(betaCI,2),length(idx),'single');
fullBetaCI(:,:,~idx) = betaCI;

 %save individual regressor sets
for iRegs = 1:length(recLabels)
    data = fullBeta(:,:,recIdx == iRegs);
    save([cPath recLabels{iRegs} 'B.mat'], 'data', 'iRegs', 'recLabels','-v7.3');
    clear data
end

%save full regressor set
save([cPath 'fullBeta.mat'], 'fullBeta', 'recIdx', 'recLabels','-v7.3');
save([cPath 'fullBetaCI.mat'], 'fullBetaCI', 'recIdx', 'recLabels','-v7.3');



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



function [grabOn,grabRel,tapOn,tapRel] = checkLevergrab(tapDur,postStimDur,grabs,release)

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
            while (grabs(Cnt+1) - release(idx)) <= 0.5 %time diff between grabs is less than 50 ms. Merge with next grab.
                
                release(idx) = []; %delete current release
                grabs(Cnt+1) = []; %delete next grab
                
                idx = find(release > grabs(Cnt),1); %next release
                if ~isempty(idx) %if there is a release that follows current grab
                    cGrab = release(idx) - grabs(Cnt); %duration of current grab
                else
                    cGrab = postStimDur - grabs(Cnt);
                end
                
                if isempty(idx) || length(grabs) <= Cnt %no more grabs/releases
                    break;
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
                grabRel(grabCnt) = idx;
            end
        end
    end
end


function dimCorr = matCorr(Vc, Vm, U)
% compute correlation between predicted and real imaging data
% Vc are real temporal components, Vm is modeled. U are spatial components

% Vm = (fullR * dimBeta)';

covVc = cov(Vc');  % S x S
covVm = cov(Vm');  % S x S
Vm = bsxfun(@minus, Vm, mean(Vm,2));
cCovV = Vm * Vc' / (size(Vc, 2) - 1);  % S x S
covP = sum((U * cCovV) .* U, 2)';  % 1 x P
varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
dimCorr = (covP ./ stdPxPy)';
