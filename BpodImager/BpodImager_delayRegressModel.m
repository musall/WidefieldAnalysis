function BpodImager_delayRegressModel(cPath,Animal,Rec)

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

%general variables
Paradigm = 'SpatialDisc';
cPath = [cPath Animal filesep Paradigm filesep Rec filesep]; %Widefield data path
sPath = ['U:\space_managed_data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
% sPath = ['/sonas-hs/churchland/hpc/home/space_managed_data/BpodImager/Animals/' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
sRate = 30; % Sampling rate of imaging in Hz
gaussShift = 3;     % inter-frame interval between regressors. Will use only every 'gaussShift' regressor and convolve with gaussian of according FHWM to reduce total number of used regressors.

preStimDur = 3;                 % Duration of trial before stimulus onset in seconds
postStimDur = 115/sRate;      % Duration of trial after stimulus onset in seconds

%other variables
mPreTime = 0.2;     % precede motor events to capture preparatory activity in seconds
mPostTime = 1;      % follow motor events for mPostStim in seconds
motorIdx = [-((mPreTime * sRate): -1 : 1) 0 (1:(mPostTime * sRate))]; %index for design matrix to cover pre- and post motor action
tapDur = 0.25;      % minimum time of lever contact, required to count as a proper grab.
noWhisk = 5;        % minimum nrs of frames before and after a lick to check for whisking. Whisks should not be counted during licking to avoid confusion.
smth = 5;           % filter length when smoothing video data
smth = smth-1+mod(smth,2); % ensure kernel length is odd
piezoLine = 2;      % channel in the analog data that contains data from piezo sensor
stimLine = 6;       % channel in the analog data that contains stimulus trigger.
ridgeFolds = 1;    %folds for cross-validation when assessing predicted variance

bhvDimCnt1 = 500;    % number of dimensions from behavioral videos that are used as regressors.
bhvDimCnt2 = bhvDimCnt1 / 5; % number of dimensions from behavioral videos that are used after 2nd SVD where videos are merged.
dims = 200; %number of dimensions from V that will be used in the model

%% load behavior data
if exist([cPath 'BehaviorVideo' filesep 'SVD_Cam 1.mat'],'file') ~= 2 %check if svd behavior exists on hdd and pull from server otherwise
    if exist([cPath 'BehaviorVideo' filesep]) ~= 2
        mkdir([cPath 'BehaviorVideo' filesep]);
    end
    copyfile([sPath 'BehaviorVideo' filesep 'SVD_Cam 1.mat'],[cPath 'BehaviorVideo' filesep 'SVD_Cam 1.mat']);
    copyfile([sPath 'BehaviorVideo' filesep 'SVD_Cam 2.mat'],[cPath 'BehaviorVideo' filesep 'SVD_Cam 2.mat']);
end

load([cPath 'BehaviorVideo' filesep 'SVD_Cam 1.mat']);

V1 = V(1:bhvDimCnt1,:)';
U1 = U;
frameTimes1 = totalFrameTimes;

load([cPath 'BehaviorVideo' filesep 'SVD_Cam 2.mat']);
V2 = V(1:bhvDimCnt1,:)';
U2 = U;
frameTimes2 = totalFrameTimes;

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
[dims, frames, trialCnt] = size(Vc);

if ~exist('bTrials','var')
    bTrials = trials;
end
bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dataset

%% find events in BPod time.
% All timestamps are relative to stimulus onset event to synchronize to imaging data later

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
    
    try
        stimTime(iTrials) = bhv.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
    catch
        stimTime(iTrials) = NaN;
    end
    
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port1In') %check for licks
        lickL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port1In;
        lickL{iTrials}(lickL{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickL{iTrials} = lickL{iTrials} - stimTime(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Port3In') %check for right licks
        lickR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Port3In;
        lickR{iTrials}(lickR{iTrials} < bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1)) = []; %dont use false licks that occured before spouts were moved in
        lickR{iTrials} = lickR{iTrials} - stimTime(iTrials);
    end
    
    leverIn(iTrials) = min(bhv.RawEvents.Trial{iTrials}.States.Reset(:)) - stimTime(iTrials); %first reset state causes lever to move in
    
    leverTimes = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
    
    stimGrab(iTrials) = leverTimes(find(leverTimes == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1) - stimTime(iTrials); %find start of lever state that triggered stimulus onset
    
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2High') %check for left grabs
        levGrabL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2High - stimTime(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1High') %check for right grabs
        levGrabR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1High - stimTime(iTrials);
    end
    
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2Low') %check for left release
        levReleaseL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2Low - stimTime(iTrials);
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1Low') %check for right release
        levReleaseR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1Low - stimTime(iTrials);
    end
    
    if ~isnan(bhv.RawEvents.Trial{iTrials}.States.Reward(1)) %check for reward state
        water(iTrials) = bhv.RawEvents.Trial{iTrials}.States.Reward(1) - stimTime(iTrials);
    end
end

 %normalize response side between -1 and 1
bhv.ResponseSide(bhv.ResponseSide == 1) = -1;
bhv.ResponseSide(bhv.ResponseSide == 2) = 1;
SessionData.ResponseSide(SessionData.ResponseSide == 1) = -1;
SessionData.ResponseSide(SessionData.ResponseSide == 2) = 1;

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

visRewardR = cell(1,trialCnt);
audRewardR = cell(1,trialCnt);

audChoiceR = cell(1,trialCnt);
visChoiceR = cell(1,trialCnt);
prevChoiceR = cell(1,trialCnt);
prevRewardR = cell(1,trialCnt);

waterR = cell(1,trialCnt);
pupilR = cell(1,trialCnt);
whiskR = cell(1,trialCnt);
faceR = cell(1,trialCnt);
bodyR = cell(1,trialCnt);
piezoR = cell(1,trialCnt);

%%
tic
for iTrials = 1:trialCnt
    %% vis/aud stim
    %regressors cover the last 2s after stimOn
    if bhv.StimType(iTrials) == 1 || bhv.StimType(iTrials) == 3 %visual or mixed stimulus
        if bhv.CorrectSide(iTrials) == 1
            lVisStimR{iTrials} = timeR(:, (preStimDur * sRate) + 1: end);
            rVisStimR{iTrials} = false(frames, postStimDur * sRate);
        else
            lVisStimR{iTrials} = false(frames, postStimDur * sRate);
            rVisStimR{iTrials} = timeR(:, (preStimDur * sRate) + 1: end);
        end
    else
        lVisStimR{iTrials} = false(frames, postStimDur * sRate);
        rVisStimR{iTrials} = false(frames, postStimDur * sRate);
    end
    
    if bhv.StimType(iTrials) == 2 || bhv.StimType(iTrials) == 3 %auditory or mixed stimulus
        if bhv.CorrectSide(iTrials) == 1
            lAudStimR{iTrials} = timeR(:, (preStimDur * sRate) + 1: end);
            rAudStimR{iTrials} = false(frames, postStimDur * sRate);
        else
            lAudStimR{iTrials} = false(frames, postStimDur * sRate);
            rAudStimR{iTrials} = timeR(:, (preStimDur * sRate) + 1: end);
        end
    else
        lAudStimR{iTrials} = false(frames, postStimDur * sRate);
        rAudStimR{iTrials} = false(frames, postStimDur * sRate);
    end
    
    if gaussShift > 1
        % subsample regressors
        lVisStimR{iTrials} = lVisStimR{iTrials}(:,1:gaussShift:end);
        rVisStimR{iTrials} = rVisStimR{iTrials}(:,1:gaussShift:end);
        lAudStimR{iTrials} = lAudStimR{iTrials}(:,1:gaussShift:end);
        rAudStimR{iTrials} = rAudStimR{iTrials}(:,1:gaussShift:end);
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
    leverInR{iTrials} = false(frames,frames);
    leverShift = round((frames/2) + leverIn(iTrials) * sRate); %timepoint in frames were lever moved in, relative to stimOnset
    
    if ~isnan(leverShift)
        if leverShift > 0 %lever moved in during the recorded trial
            leverInR{iTrials}(:, 1: frames - leverShift) = timeR(:, leverShift+1:end);
        else %lever was present before data was recorded
            leverInR{iTrials}(:, abs(leverShift) + 1 : end) = timeR(:,  1: frames + leverShift);
        end
    end
    
    if gaussShift > 1
        % subsample regressors
        leverInR{iTrials} = leverInR{iTrials}(:,1:gaussShift:end);
    end     

    %% choice and reward
    visRewardR{iTrials} = logical((timeR .* (bhv.StimType(iTrials)  == 1)) * bhv.Rewarded(iTrials)); %visual trial, 1 for reward
    audRewardR{iTrials} = logical((timeR .* (bhv.StimType(iTrials)  == 2)) * bhv.Rewarded(iTrials)); %audio trial, 1 for reward
    
    visChoiceR{iTrials} = single(timeR .* bhv.ResponseSide(iTrials) .* (bhv.StimType(iTrials)  == 1)); %visual choice regressor. Takes -1 for left and 1 for right responses
    audChoiceR{iTrials} = single(timeR .* bhv.ResponseSide(iTrials) .* (bhv.StimType(iTrials)  == 2)); %auditory choice regressor. Takes -1 for left and 1 for right responses
    
    if iTrials == 1 %don't use first trial
        prevRewardR{iTrials} = NaN(size(timeR));
        prevChoiceR{iTrials} = NaN(size(timeR));
        
    else %for all subsequent trials
        prevChoiceR{iTrials} = single(timeR .* SessionData.ResponseSide(bTrials(iTrials)-1)); % check if previous trial was rewarded (1), punished (-1). Goes to 0 if not performed
        if ~SessionData.Rewarded(bTrials(iTrials)-1) %last trial was not rewarded
            prevChoiceR{iTrials} = zeros(size(timeR),'single');
            prevRewardR{iTrials} = zeros(size(timeR),'single');
        else
            prevRewardR{iTrials} = single(timeR);
        end
    end
    
    if gaussShift > 1
        % subsample regressors
        visRewardR{iTrials} = visRewardR{iTrials}(:,1:gaussShift:end);
        audRewardR{iTrials} = audRewardR{iTrials}(:,1:gaussShift:end);
        prevRewardR{iTrials} = prevRewardR{iTrials}(:,1:gaussShift:end);
        visChoiceR{iTrials} = visChoiceR{iTrials}(:,1:gaussShift:end);
        audChoiceR{iTrials} = audChoiceR{iTrials}(:,1:gaussShift:end);
        prevChoiceR{iTrials} = prevChoiceR{iTrials}(:,1:gaussShift:end);
    end
    
    %determine timepoint of reward given
    waterR{iTrials} = false(frames, sRate);
    if ~isnan(water(iTrials)) && ~isempty(water(iTrials))
        waterOn = round((frames/2) + water(iTrials) * sRate); %timepoint in frames when reward was given
        waterR{iTrials}(:, 1: size(timeR,2) - waterOn) = timeR(:, waterOn+1:end);
    end
    
    if gaussShift > 1
        % subsample regressors
        waterR{iTrials} = waterR{iTrials}(:,1:gaussShift:end);
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
    
    %% pupil/video regressors
    frameFile = dir([cPath 'BehaviorVideo' filesep '*frameTimes_'  num2str(bTrials(iTrials),'%04i') '*_1.mat']);
    if exist([cPath 'BehaviorVideo' filesep 'faceVars_' int2str(bTrials(iTrials)) '.mat'],'file') ~= 2 || isempty(frameFile)  %check if files exists on hdd and pull from server otherwise
        if exist([cPath 'BehaviorVideo' filesep]) ~= 2
            mkdir([cPath 'BehaviorVideo' filesep]);
        end
        frameFile = dir([sPath 'BehaviorVideo' filesep '*frameTimes_'  num2str(bTrials(iTrials),'%04i') '*_1.mat']);
        
        %get results from pupil analysis
        copyfile([sPath filesep 'BehaviorVideo' filesep 'faceVars_' int2str(bTrials(iTrials)) '.mat'],[cPath filesep 'BehaviorVideo' filesep 'faceVars_' int2str(bTrials(iTrials)) '.mat']);
        % copy frametimes from server to local
        copyfile([sPath filesep 'BehaviorVideo' filesep frameFile.name],[cPath filesep 'BehaviorVideo' filesep frameFile.name]);
    end
    
    load([cPath 'BehaviorVideo' filesep 'faceVars_' int2str(bTrials(iTrials)) '.mat'])
    load([cPath 'BehaviorVideo' filesep frameFile.name])
    bhvStart(iTrials) = frameTimes(1) * 86400;
    
    %absolute trial onset time for bpod - use to synchronize video to behavioral events
    trialOn = bhv.TrialStartTime(iTrials) + (stimTime(iTrials) - preStimDur);
    frameTimes = (frameTimes  * 86400) - trialOn;
    trialOn = find(frameTimes > 0,1); %trial on in frames
    bhvFrameRate = round(1/median(diff(frameTimes)));
    trialOn = round(trialOn * (sRate / bhvFrameRate)); %trial on in resampled frames
    
    %resample to match imaging data framerate
    try
        pupil = mean(eyeVars.axes,2); %pupil diameter
        idx = zscore(diff(pupil)) > 0.5; %index for blinks
        idx(1) = false; idx(end) = false; %first and last datapoints should not be NaN
        t = 1:length(pupil);
        pupil(idx) = interp1(t(~idx), pupil(~idx), t(idx)); % do interpolation to bridge gaps in pupil trace
        pupil = [repmat(pupil(1),21,1); pupil; repmat(pupil(end),21,1)]; %add some padding on both sides to avoid edge effects when resampling
        pupil = resample(pupil, sRate, bhvFrameRate); %resample to match imaging data framerate
        offset = ceil(21 * sRate / bhvFrameRate); %required offset to remove padding after resampling
        pupil = pupil(offset : end - offset); %remove padds
        pupil = smooth(pupil,'rlowess'); %do some smoothing
        idx = trialOn:trialOn + (frames - 1); %index for data that matches the widefield
        pupilR{iTrials} = single(pupil(idx)); %only use trial-relevant frames
    catch
        pupilR{iTrials} = NaN(frames,1,'single'); %can't use this trial
    end
        
    % get snoutmotion, using the same method as before
    try
        snoutMotion = [repmat(snoutMotion(1),21,1); snoutMotion'; repmat(snoutMotion(end),21,1)]; %add some padding on both sides to avoid edge effects when resampling
        snoutMotion = resample(snoutMotion, sRate, bhvFrameRate); %resample to match imaging data framerate
        offset = ceil(21 * sRate / bhvFrameRate); %required offset to remove padding after resampling
        snoutMotion = snoutMotion(offset : end - offset); %remove padds
        snoutMotion = smooth(snoutMotion,'rlowess'); %do some smoothing
        idx = trialOn:trialOn + (frames - 1); %index for data that matches the widefield
        whiskR{iTrials} = single(snoutMotion(idx)); %only use trial-relevant frames
    catch
        whiskR{iTrials} = NaN(frames,1,'single'); %can't use this trial
    end
    
    %get svd data from cam1 - this is usually the one that sees the animals face
    trialOn = bhv.TrialStartTime(iTrials) + (stimTime(iTrials) - preStimDur);
    cTimes = frameTimes1 * 86400 - trialOn;
    trialOn = find(cTimes > 0,1); %trial on in frames
    bhvFrameRate = round(1/median(diff(cTimes)));
    idx = trialOn : trialOn + ceil(frames*(bhvFrameRate / sRate)); %index for data that matches the widefield after resampling
    try
        faceMotion = V1(idx, 1:bhvDimCnt1);
        faceMotion = [repmat(faceMotion(1,:),21,1); faceMotion; repmat(faceMotion(end,:),21,1)]; %add some padding on both sides to avoid edge effects when resampling
        faceMotion = resample(double(faceMotion), sRate, bhvFrameRate); %resample to match imaging data framerate
        faceMotion = conv2(faceMotion',ones(1,smth)/smth,'same')'; %smooth trace with moving average of 'smth' points
        faceMotion = faceMotion(offset : end - offset,:); %remove padds
        faceR{iTrials} = single(faceMotion(1:frames,:)); %only use trial-relevant frames
    catch
        faceR{iTrials} = NaN(frames,bhvDimCnt1, 'single'); %can't use this trial
    end
    
    %get svd data from cam2 - this is usually the one that sees the animal from the bottom
    trialOn = bhv.TrialStartTime(iTrials) + (stimTime(iTrials) - preStimDur);
    cTimes = frameTimes2 * 86400 - trialOn;
    trialOn = find(cTimes > 0,1); %trial on in frames
    bhvFrameRate = round(1/median(diff(cTimes)));
    idx = trialOn : trialOn + ceil(frames*(bhvFrameRate / sRate)); %index for data that matches the widefield after resampling
    try
        bodyMotion = V2(idx, 1:bhvDimCnt1);    
        bodyMotion = [repmat(bodyMotion(1,:),21,1); bodyMotion; repmat(bodyMotion(end,:),21,1)]; %add some padding on both sides to avoid edge effects when resampling
        bodyMotion = resample(double(bodyMotion), sRate, bhvFrameRate); %resample to match imaging data framerate
        bodyMotion = conv2(bodyMotion',ones(1,smth)/smth,'same')'; %smooth trace with moving average of 'smth' points
        bodyMotion = bodyMotion(offset : end - offset,:); %remove padds
        bodyR{iTrials} = single(bodyMotion(1:frames,:)); %only use trial-relevant frames
    catch
        bodyR{iTrials} = NaN(frames, bhvDimCnt1, 'single'); %can't use this trial
    end
    
    %% piezo sensor information
    if exist([cPath 'Analog_'  num2str(trials(iTrials)) '.dat'],'file') ~= 2  %check if files exists on hdd and pull from server otherwise
        copyfile([sPath 'Analog_'  num2str(trials(iTrials)) '.dat'],[cPath 'Analog_'  num2str(trials(iTrials)) '.dat']);
    end
    
    [~,Analog] = Widefield_LoadData([cPath 'Analog_'  num2str(trials(iTrials)) '.dat'],'Analog'); %load analog data    
    stimOn = find(diff(double(Analog(stimLine,:)) > 1500) == 1); %find stimulus onset in current trial
    Analog(1,round(stimOn + (postStimDur * 1000) - 1)) = 0; %make sure there are enough datapoints in analog signal
    piezoR{iTrials} = Analog(piezoLine,round(stimOn - (preStimDur * 1000)) : round(stimOn + (postStimDur * 1000) - 1)); % data from piezo sensor. Should encode animals hindlimb motion.
        
    % give some feedback over progress
    if rem(iTrials,50) == 0
        fprintf(1, 'Current trial is %d out of %d\n', iTrials,trialCnt);
        toc
    end
end

%% rebuild piezo sensor to get proper design matrix
wTrace = double(cat(1,whiskR{:}));
idx = ~isnan(wTrace);
wTrace(idx) = smooth(wTrace(idx));
wTrace = (wTrace - prctile(wTrace,1))./ nanstd(wTrace); %minimum values are at 0, signal in standard deviation units
wTrace = wTrace > 2; %take activity above 4 SDUs as indicator for whisking
wTrace = diff([0; wTrace]) == 1; %find event onsets
wTrace = reshape(wTrace,[],trialCnt);

for iTrials = 1:trialCnt
    
    trace = logical(histcounts(find(wTrace(:,iTrials)),0:bhvFrameRate/sRate:(bhvFrameRate/sRate)*frames))'; %resample to imaging frame rate. This is the zero lag regressor.    
    licks = round(([lickL{iTrials} lickR{iTrials}] + preStimDur) * sRate);  %lick times in frames. Dont check for whisking there because licking may move the snout around as well.
    if ~isempty(licks)
        lickZone = reshape(bsxfun(@plus, licks', -noWhisk : noWhisk),1,[]); %done use whisks that are too close to licks.
        trace(lickZone) = false;
    end
    
    % create full design matrix
    cIdx = bsxfun(@plus,find(trace),motorIdx);
    cIdx(cIdx < 1) = 0;
    cIdx(cIdx > frames) = frames;
    cIdx = bsxfun(@plus,cIdx,(0:frames:frames*length(motorIdx)-1));
    cIdx(cIdx < 1) = frames;
    cIdx(cIdx > (frames * length(motorIdx))) = frames * length(motorIdx);
    
    whiskR{iTrials} = false(frames, length(motorIdx));
    whiskR{iTrials}(cIdx(:)) = true;
    whiskR{iTrials}(end,:) = false; %don't use last timepoint of design matrix to avoid confusion with indexing.
    
    if gaussShift > 1
        whiskR{iTrials} = whiskR{iTrials}(:,1:gaussShift:end);
    end
end
clear wTrace idx cIdx

%% rebuild piezo sensor to get proper design matrix
pTrace = double(cat(2,piezoR{:}));
pTrace = smooth(zscore(pTrace),100);
pTrace = pTrace > 0.5; %threshold normalized sensor data
pTrace = diff([0;pTrace]) == 1; %find event onsets
pTrace = reshape(pTrace,[],trialCnt);

for iTrials = 1:trialCnt
    trace = logical(histcounts(find(pTrace(:,iTrials)), 0:1000/sRate:(1000/sRate)*frames))'; %resample to imaging frame rate. This is the zero lag regressor.
    cIdx = bsxfun(@plus,find(trace),motorIdx);
    cIdx(cIdx < 1) = 0;
    cIdx(cIdx > frames) = frames;
    cIdx = bsxfun(@plus,cIdx,(0:frames:frames*length(motorIdx)-1));
    cIdx(cIdx < 1) = frames;
    cIdx(cIdx > (frames * length(motorIdx))) = frames * length(motorIdx);
    
    piezoR{iTrials} = false(frames, length(motorIdx));
    piezoR{iTrials}(cIdx(:)) = true;
    piezoR{iTrials}(end,:) = false; %don't use last timepoint of design matrix to avoid confusion with indexing.
    
    if gaussShift > 1
        piezoR{iTrials} = piezoR{iTrials}(:,1:gaussShift:end);
    end
end
clear pTrace cIdx

%% reshape regressors, make design matrix and indices for regressors that are used for the model
timeR = logical(diag(ones(1,frames)));
timeR = timeR(:, 1:preStimDur*sRate); %only use baseline part of timeR to avoid redundancy with stimulus regressors

if gaussShift > 1
    timeR = timeR(:,1:gaussShift:end); %subsample regressor matrix
end
timeR = repmat(timeR,length(trials),1);

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
prevRewardR = cat(1,prevRewardR{:});

visChoiceR = cat(1,visChoiceR{:});
audChoiceR = cat(1,audChoiceR{:});
prevChoiceR = cat(1,prevChoiceR{:});

waterR = cat(1,waterR{:});

pupilR = cat(1,pupilR{:});
camIdx = ~isnan(pupilR(:,1)); %index for trials where bhv frames were available
pupilR(camIdx,:) = zscore(pupilR(camIdx,:));

piezoR = cat(1,piezoR{:});
whiskR = cat(1,whiskR{:});

%% combine body and face data and use svd to merge into one set of predictors (to reduce redundancy between video regressors)
videoR = [cat(1,faceR{:}) cat(1,bodyR{:})]; clear faceR bodyR
videoR = videoR(camIdx,:);

[vidV, mVidV, stdVidV, vidU] = compressVideoRegs(videoR,bhvDimCnt2); %compress into one set of regress and zscore so it can be used in GLM
vidR = NaN(length(camIdx),bhvDimCnt2,'single');
vidR(camIdx,:) = vidV; clear vidV
save([cPath filesep 'vidR.mat'], 'vidR', 'mVidV' ,'stdVidV', 'vidU');

%% create full design matrix
fullR = [timeR lGrabR lGrabRelR rGrabR rGrabRelR lLickR rLickR leverInR lVisStimR rVisStimR lAudStimR rAudStimR visRewardR audRewardR ...
    prevRewardR (visChoiceR+audChoiceR) prevChoiceR waterR piezoR whiskR pupilR vidR];

trialIdx = isnan(mean(fullR,2)); %don't use first trial or trials that failed to contain behavioral video data
fullR(trialIdx,:) = []; %clear bad trials
fprintf(1, 'Rejected %d/%d trials for NaN regressors\n', sum(trialIdx) / frames,trialCnt);

idx = nansum(abs(fullR)) < 10; %reject regressors that are too sparse
fullR(:,idx) = []; %clear empty regressors
fprintf(1, 'Rejected %d/%d empty regressors\n', sum(idx),length(idx));

% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
recLabels = {
    'time' 'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick' 'leverIn' 'lVisStim' 'rVisStim' ...
    'lAudStim' 'rAudStim' 'visReward' 'audReward' 'prevReward' 'visChoice' 'audChoice' 'prevChoice' ...
    'water' 'piezo' 'whisk' 'pupil' 'bhvVideo'};

%index to reconstruct different response kernels
recIdx = [
    ones(1,size(timeR,2))*find(ismember(recLabels,'time')) ... 
    ones(1,size(lGrabR,2))*find(ismember(recLabels,'lGrab')) ... 
    ones(1,size(lGrabRelR,2))*find(ismember(recLabels,'lGrabRel')) ... 
    ones(1,size(rGrabR,2))*find(ismember(recLabels,'rGrab')) ...
    ones(1,size(rGrabRelR,2))*find(ismember(recLabels,'rGrabRel')) ...
    ones(1,size(lLickR,2))*find(ismember(recLabels,'lLick')) ...
    ones(1,size(rLickR,2))*find(ismember(recLabels,'rLick')) ...
    ones(1,size(leverInR,2))*find(ismember(recLabels,'leverIn')) ...
    ones(1,size(lVisStimR,2))*find(ismember(recLabels,'lVisStim')) ...
    ones(1,size(rVisStimR,2))*find(ismember(recLabels,'rVisStim')) ...
    ones(1,size(lAudStimR,2))*find(ismember(recLabels,'lAudStim')) ...
    ones(1,size(rAudStimR,2))*find(ismember(recLabels,'rAudStim')) ...
    ones(1,size(visRewardR,2))*find(ismember(recLabels,'visReward')) ...
    ones(1,size(audRewardR,2))*find(ismember(recLabels,'audReward')) ...
    ones(1,size(prevRewardR,2))*find(ismember(recLabels,'prevReward')) ...
    ones(1,size(visChoiceR,2))*find(ismember(recLabels,'visChoice')) ...
    ones(1,size(audChoiceR,2))*find(ismember(recLabels,'audChoice')) ...
    ones(1,size(prevChoiceR,2))*find(ismember(recLabels,'prevChoice')) ...
    ones(1,size(waterR,2))*find(ismember(recLabels,'water')) ...
    ones(1,size(piezoR,2))*find(ismember(recLabels,'piezo')) ...
    ones(1,size(whiskR,2))*find(ismember(recLabels,'whisk')) ...
    ones(1,size(pupilR,2))*find(ismember(recLabels,'pupil')) ...
    ones(1,size(vidR,2))*find(ismember(recLabels,'bhvVideo'))];

%% apply gaussian filter to design matrix if using sub-sampling and check if it is sufficently orthogonal
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

[Q, R] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize design matrix

pupilIdx = recIdx(~idx) == find(ismember(recLabels,'pupil'));
fullR(:,pupilIdx) = Q(:,pupilIdx); %orthogonalize pupil regressors

vidIdx = recIdx(~idx) == find(ismember(recLabels,'bhvVideo'));
fullR(:,vidIdx) = Q(:,vidIdx); %orthogonalize video regressors

if ~fullColRank(R) %check if design matrix is full rank
    error('Design matrix is rank-defficient')
end

% save some results
save([cPath filesep 'regData.mat'], 'fullR', 'idx' ,'trialIdx', 'recIdx', 'recLabels');

%clear individual regressors
clear stimR lGrabR lGrabRelR rGrabR rGrabRelR waterR lLickR rLickR leverInR ...
    lVisStimR rVisStimR lAudStimR rAudStimR visRewardR audRewardR prevRewardR visChoiceR audChoiceR ...
    prevChoiceR pupilR bhvVideo piezoR whiskR

%% run ridge regression in low-D and find error term
U = arrayShrink(U,mask);
Vc = reshape(Vc,dims,[]);
Vc(:,trialIdx) = [];
Vc = bsxfun(@minus, Vc, mean(Vc, 2)); %make sure Vc is zero-mean

[ridgeVals, dimBeta] = ridgeMML(Vc', fullR); %get beta weights
save([cPath 'dimBeta.mat'], 'dimBeta', 'ridgeVals')
fprintf('Mean ridge penalty: %f\n', mean(ridgeVals));

%% collect beta maps and reconstruct kernels
beta = dimBeta * U'; %compute beta weight maps
beta = arrayShrink(gather(beta)',mask,'split'); %spatially reconstruct beta maps
fullBeta = NaN(size(beta,1),size(beta,2),length(idx),'single');
fullBeta(:,:,~idx) = beta; %full matrix to include regressors that were rejected before
save([cPath 'fullBeta.mat'], 'fullBeta', 'recIdx', 'recLabels', 'motorIdx', 'gaussShift','-v7.3'); %save full regressor set

 %save individual regressor sets
for iRegs = 1:length(recLabels)
    
    data = fullBeta(:,:,recIdx == iRegs);
    data = reshape(data,[],size(data,3));
    
    temp = 1:gaussShift:sum(recIdx == iRegs)*gaussShift;
    temp = zeros(temp(end),numel(mask));
    temp(1:gaussShift:end,:) = data';
    temp = smoothCol(temp,gaussShift*2,'gauss')';
    data = single(reshape(temp,size(mask,1),size(mask,2),[]));
        
    save([cPath recLabels{iRegs} 'B.mat'], 'data', 'iRegs','recLabels','gaussShift','-v7.3');
    clear data temp
end

%% test predictive power of individual regressors through cross-validation
disp(size(Vc));
predTrial = zeros(length(unique(recIdx))+1,frames);

Cnt = 0;
for iRegs = 0 : length(recLabels)
    
    Cnt = Cnt+1;
    fprintf('Current regressor is %d out of %d\n', Cnt,length(unique(recIdx))+1);
    
    cIdx = recIdx(~idx) == iRegs; %index for reduced model.
    fakeR = fullR; %copy  design matrix to shuffle current regressor set
    for iCol = find(cIdx)
        fakeR(:,iCol) = fullR(randperm(size(fullR,1)),iCol); %shuffle current regressors in time
    end
%     fakeR(:,cIdx) = [];
    
    [dimRsq, meanRsq, trialRsq, cMap, cMovie] = Widefield_crossValModel(Vc,fakeR,U,ridgeFolds,frames,ispc);
       
    if iRegs == 0 %first regressor, collect correlation map
        cMovie = sqrt(cMovie.^2);
        fullRsq = median(dimRsq);
        fullMovie = cMovie;
    else
        cMovie = fullMovie - sqrt(cMovie.^2); %any other regressor: subtract shuffled from full model to identify areas of reduced correlation
    end
    pMovie = reshape(corrTest(cMovie(:),single(trialCnt)),[],frames);
    [rankP(Cnt), rankH(Cnt)] = signrank(dimRsq, fullRsq, 'tail', 'left'); %check for significant reduction vs full model
    
    if iRegs == 0
        save([cPath 'fullcorr.mat'], 'meanRsq', 'dimRsq', 'trialRsq', 'cMap', 'cMovie', 'pMovie', 'iRegs','recLabels','gaussShift','-v7.3');
    else
        save([cPath recLabels{iRegs} 'corr.mat'], 'meanRsq', 'dimRsq', 'trialRsq', 'cMap', 'cMovie', 'pMovie', 'iRegs','recLabels','gaussShift','-v7.3');
    end
    fprintf('Finished. Current reg: %d ; pVal is %d\n', Cnt,rankP(Cnt));

end
save([cPath recLabels{iRegs} 'regRankTest.mat'], 'rankP', 'rankH');


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

function [vidV, meanV, stdV, vidU] = compressVideoRegs(vidR,dimCnt)
% perform svd to further compres video regressors and reduce redundancy
% between cameras

[vidV,S,vidU]  = svd(vidR, 'econ');
vidU = vidU(:,1:dimCnt)'; %the resulting U is nSVD x regs 
vidV = vidV*S; %convolve V and S to stick to the usual analysis. This is frames x nSVD and is used as regressor in the GLM.
vidV = vidV(:,1:dimCnt); %get request nr of regressors

%zscore video regressor, return mean and std for later reconstruction
meanV = mean(vidV);
stdV = std(vidV);
vidV = bsxfun(@minus,vidV,meanV);
vidV = bsxfun(@rdivide,vidV,stdV);