function timeDiff = Widefield_checkMissedTrials(cPath)
% this code was an attempt to align imaging and behavioral trials when barcodes were not correctly recorded and trials are missing.
% doesnt work very well so far - would need more work to be trusted.


if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

%% check if .dat or .mj2 files are present. Use slower splitvideo code for latter case.
rawCheck = dir([cPath filesep 'Frames*dat']);
vidCheck = dir([cPath filesep 'frameTimes*mat']);
analogCheck = dir([cPath filesep 'Analog_*.dat']);

if size(rawCheck,1) == size(analogCheck,1)
    loadRaw = true;
elseif size(vidCheck,1) == size(analogCheck,1)
    loadRaw = false;
else
    error('Unequal number of imaging and analog data files. Aborted')
end

%%

if loadRaw
else
    for iFiles = 1:length(vidCheck)
        
        load([cPath vidCheck(iFiles).name])
        frameTimes = frameTimes * (86400); %extract frame times from header and convert to millisecond timestamps
        trialOn(iFiles) = frameTimes(1);
        trialOff(iFiles) = frameTimes(end);
        
    end
end

[~, b] = sort(trialOn);
trialOn = trialOn(b);
trialOff = trialOff(b);
trials = 1:length(trialOn);

%%
figure
timeDiff = ((trialOn(2:end) - trialOff(1:end-1)));
histogram(timeDiff); axis square; xlabel('time (s)'); ylabel('# trials');

%% get behavioral data to get stimulus onset times
bhvFile = strsplit(cPath,filesep);
bhvFile = dir([cPath bhvFile{end-3}(1:5) '*' bhvFile{end-2} '*.mat']);
load([cPath bhvFile.name]);
bhvOn = SessionData.TrialStartTime .* (86400);
bTrials = 1:length(bhvOn);

%% 
% for iRuns = 1 : length(bhvOn)-length(trialOn)
%     [~,a] = max(abs(diff(trialOn - bhvOn(1:length(trialOn)))))
%     bhvOn(a) = [];
% end



%%
Cnt = 0;
for iTrials = 1:length(trialOn)
    while abs(trialOn(iTrials) - bhvOn(iTrials)) > prctile(timeDiff,98)
        bhvOn(iTrials) = [];
        bTrials(iTrials) = [];
        Cnt = Cnt +1;
    end
end


end