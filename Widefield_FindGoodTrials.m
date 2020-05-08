function [trials, frameCnt, nTrials, bTrials] = Widefield_FindGoodTrials(opts)
% Code to identify trials that should be used for imaging analysis. 
% Reject trials with incomplete frame count and add indicator for touch
% events during baseline to behavioral data.

%% find all trial files and check barodes
recs = dir([opts.fPath filesep opts.fName '*']);
trials = zeros(1,length(recs));
for iRecs = 1:length(recs)
    a = textscan(recs(iRecs).name,'%s%f%s','Delimiter','_');
    trials(iRecs) = a{2};
end
trials = sort(trials);
       
%% load analog data to check trial barcodes. Remove unassigned trials from imaging data.
bTrials = zeros(1,length(trials));
frameCnt = zeros(1,length(trials));
stimCnt = zeros(1,length(trials));

for iTrials = 1:length(trials)
    cFile = [opts.fPath filesep 'Analog_' num2str(trials(iTrials)) '.dat']; %analog file to be read
    [~,Analog] = Widefield_LoadData(cFile,'Analog'); %load analog data
    
    if opts.loadRaw
        cFile = [opts.fPath filesep opts.fName '_' num2str(trials(iTrials)) '.dat']; %imaging file to be read
        header = Widefield_LoadData(cFile,'Analog');
        frameCnt(iTrials) = header(end);
    else
        cFile = [opts.fPath filesep 'frameTimes_' num2str(trials(iTrials)) '.mat']; %file with timestamps
        temp = load(cFile);
        frameCnt(iTrials) = length(temp.frameTimes);
    end
    
    if opts.barcodeCheck
        try
            bTrials(iTrials) = Widefield_readBarcodes(double(Analog(opts.barcodeLine,:))./4095*5, 3, 6); %Bpod TrialNr, encoded from barcode signal
        catch
            bTrials(iTrials) = -1; % don't use this trial if barcode can't be recovered
        end
    else
        bTrials(iTrials) = iTrials;  %skip this control if barcodes are not readable for some reason.
    end
    
    % check for missing stimulus triggers
    stimOn = diff(double(Analog(opts.stimLine,:)) > 1500) == 1;
    if find(stimOn) < length(stimOn) - 250 %needs at least 250ms of analog data after stimOnset
        stimCnt(iTrials) = sum(stimOn);
    end    
end
clear Analog

if isfield(opts,'minFrameCnt')
    minTrials = opts.minFrameCnt;
else
    minTrials = prctile(frameCnt,90)*.75;
end

% only use trials that have a correct stimulus trigger
if any(sum(stimCnt(frameCnt >= minTrials) == 0))
    fprintf('Rejected %d/%d trials for missing trigger signal \n', sum(stimCnt(frameCnt >= minTrials) == 0), length(bTrials))
end

% only use trials that have a complete framecount and working barcode
if any(frameCnt < round(minTrials))
    fprintf('Rejected %d/%d trials for insufficient acquired frames \n', sum(frameCnt < round(minTrials)), length(bTrials))
end

ind = (frameCnt < minTrials) | bTrials <= 0 | stimCnt ~= 1; 
trials(ind) = [];
bTrials(ind) = [];
frameCnt(ind) = [];

%% load behavior data, check stim delay and add indicies for potentially unfit baseline activity
cFile = dir([opts.fPath filesep opts.animal '_' opts.paradigm '*.mat']);
cFile = strtrim(cFile.name);
load([opts.fPath filesep cFile]); %load behavior data

ind = bTrials > SessionData.nTrials; %check for imaging trials that are not in the bhv file
trials(ind) = [];
bTrials(ind) = [];
frameCnt(ind) = [];
if sum(ind) > 0
    fprintf('Rejected %d/%d trials that are not present in Bpod data \n', sum(ind), length(ind))
end

badTrial = false(1,SessionData.nTrials);
SessionData.TouchCnt = zeros(1,length(SessionData.Rewarded));
for iTrials = bTrials
    
    stimDelay = [];
    stimTime = SessionData.RawEvents.Trial{iTrials}.States.RunningStimulus;
    % check for delayed stimulus responses
    if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Wire3High')
        stimDelay = SessionData.RawEvents.Trial{iTrials}.Events.Wire3High - stimTime(1); %recorded stimOn minus presumed stimOn
    end
    
    stimDelay(stimDelay < 0) = [];
    if any(stimDelay > 0.1) || isempty(stimDelay) % %delay above 100ms is not ok
        badTrial(iTrials) = true;
    end
    
    leverTouch = [];
    if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Wire1Low')
        leverTouch = SessionData.RawEvents.Trial{iTrials}.Events.Wire1Low;
    end
    if isfield(SessionData.RawEvents.Trial{iTrials}.Events,'Wire2Low')
        leverTouch = [leverTouch SessionData.RawEvents.Trial{iTrials}.Events.Wire2Low];
    end
    SessionData.TouchCnt(iTrials) = sum(((sort(leverTouch) - stimTime(1))) < 0 &  ((sort(leverTouch) - stimTime(1))) >= -(SessionData.Settings.preStimDelay + 1));
end

if any(badTrial)
    fprintf('Found %d/%d bad trials due to stimOn delays above 100ms.\n',sum(badTrial),length(badTrial))
    frameCnt(ismember(bTrials,find(badTrial))) = [];
    trials(ismember(bTrials,find(badTrial))) = [];
    bTrials(ismember(bTrials,find(badTrial))) = [];
end

save([opts.fPath filesep cFile], 'SessionData'); %save behavior data
nTrials = SessionData.nTrials;