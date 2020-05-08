function BpodImager_regressModel(Animal,Rec)


Paradigm = 'SpatialDisc';
cPath = ['H:\BpodImager\Animals\' Animal filesep Paradigm filesep Rec]; %Widefield data path
sPath = ['U:\space_managed_data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec ]; %server data path
sRate = 40;         % Sampling rate of imaging in Hz
preStimDur = 2;    % Duration of trial before stimulus onset in seconds
postStimDur = 2;    % Duration of trial after stimulus onset in seconds
dimCnt = 250;       % number of dimensions used for analysis
gaussWidth = 200;
mPretime = 1;     % precede self-initiated events in case a neuron predicts animal motor action
sPretime = 0;

%% load data
load([cPath '\Vc.mat'])
load([cPath '\mask.mat'])
load([cPath '\' Animal '_opts'])
load([cPath filesep ls([cPath filesep Animal '_' Paradigm '*.mat'])]); %load behavior data
bhv = selectBehaviorTrials(SessionData,trials); %only use completed trials that are in the Vc dataset

%% get time point of stimulus onset from analog data
stimOn = BpodImager_recoverStimOn(cPath,sPath,trials);

%% find events in BPod time. 
% All timestamps are relative to stimulus onset event to synchronize to imaging data later

% pre-allocate vectors
lickL = cell(1,length(trials));
lickR = cell(1,length(trials));
leverIn = NaN(1,length(trials));
spoutBias = NaN(1,length(trials));
visStim = NaN(1,length(trials));
audStim = NaN(1,length(trials));
levGrabL = cell(1,length(trials));
levGrabR = cell(1,length(trials));
levReleaseL = cell(1,length(trials));
levReleaseR = cell(1,length(trials));

for iTrials = 1:length(trials)
    
    if bhv.StimType(iTrials) == 1 || bhv.StimType(iTrials) == 3 %visual or mixed stimulus
        visStim(iTrials) = bhv.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
    end
    if bhv.StimType(iTrials) == 2 || bhv.StimType(iTrials) == 3 %auditory or mixed stimulus
        audStim(iTrials) = bhv.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
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
    spoutBias(iTrials) = bhv.TrialSettings(iTrials).ServoPos(2) - bhv.TrialSettings(iTrials).ServoPos(1); %negative values are left bias, positive are right bias
       
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1High') %check for left grabs
        levGrabL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1High - stimTime;
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2High') %check for right grabs
        levGrabR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2High - stimTime;
    end
        
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire1Low') %check for left release
        levReleaseL{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire1Low - stimTime;
    end
    if isfield(bhv.RawEvents.Trial{iTrials}.Events,'Wire2Low') %check for right release
        levReleaseR{iTrials} = bhv.RawEvents.Trial{iTrials}.Events.Wire2Low - stimTime;
    end
end

%% build regressors - create design matrix based on event times
% spoutMove = round(bhv.RawEvents.Trial{iTrials}.States.MoveSpout(1) - bhv.RawEvents.Trial{iTrials}.Events.Wire3High(1),3); %round to ms
spoutMove = 1; %spouts move in 1s after stimOn
lickWindow = postStimDur - (spoutMove - mPretime); %time window to be covered by lick regressors

tBin = round((gaussWidth/2) / (1/sRate*1000));
sigma = gaussWidth / (2*sqrt(2*log(2))); %compute std based on required full width half max
x = 0 : (1/sRate)*1000 : 3*sigma; %compute full width of the gaussian
kernel = normpdf(sort(unique([-x x])),0,sigma); %gaussian kernel
kernel = kernel / max(kernel);

%create basic time regressors and convolve with gaussian
timeR = zeros(size(Vc,2),size(Vc,2)/tBin);
idx = sub2ind(size(timeR), 1:tBin:size(timeR,2)*tBin, 1:size(timeR,2));
timeR(idx) = 1;
timeR = conv2(timeR',kernel,'same')'; %convolve with gaussian kernel

for iTrials = 1:length(trials)
    %% vis/aud stim
    %regressors cover the last 2s after stimOn
    vStimR = zeros(size(Vc,2),(postStimDur * sRate) / tBin);
    aStimR = zeros(size(Vc,2),(postStimDur * sRate) / tBin);
    
    idx = (1 : tBin : size(vStimR,2)*tBin) + (preStimDur * sRate);
    idx = sub2ind(size(vStimR), idx, 1:size(vStimR,2));
    
    if ~isnan(visStim(iTrials)) && ~isempty(visStim(iTrials))
        vStimR(idx) = 1;
        vStimR = conv2(vStimR',kernel,'same')';
    end
    if ~isnan(audStim(iTrials)) && ~isempty(audStim(iTrials))
        aStimR(idx) = 1;
        aStimR = conv2(aStimR',kernel,'same')';
    end

    %% lick regressors
    lLickR = zeros(size(Vc,2),(lickWindow * sRate) / tBin);
    rLickR = zeros(size(Vc,2),(lickWindow * sRate) / tBin);
    
    for iRegs = 0:size(lLickR,2)-1
        
        licks = lickL{iTrials} - (mPretime - (iRegs * tBin * 1/sRate));
        lLickR(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
        
        licks = lickR{iTrials} - (mPretime - (iRegs * tBin * 1/sRate));
        rLickR(logical(histcounts(licks,-preStimDur:1/sRate:postStimDur)),iRegs+1) = 1;
        
    end
    
    lLickR = conv2(lLickR',kernel,'same')'; %convolve with gaussian kernel
    rLickR = conv2(rLickR',kernel,'same')'; %convolve with gaussian kernel

    %% lever in
    leverShift = round((leverIn(iTrials) + preStimDur) * sRate); %timepoint in frames were lever moved in, relative to stimOnset
    timeR
    
    
%% spout bias
spoutBias = spoutBias / max(abs(spoutBias));
sBiasR = timeR .* spoutBias(iTrials);


%%

end

