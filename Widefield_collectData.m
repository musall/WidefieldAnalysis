function [mov, frameCnt, blueFrameTimes, hemoFrameTimes, stimOn, blueAvg, hemoAvg, opts] = Widefield_collectData(cPath, Animal, Rec)
% code to collect all data from a given imaging experiment in a larger
% matrix 'mov' that can be used to dimensionality reduction code.
% mov is pixels x frames x channels and normalized by the average over all
% used frames.
% channels are either blue (1) or violet (2) which contains intrinsic signals.
% frameCnt is a vector that contains the number of used frames per trial
% blue/hemoFrameTimes are relative timestamps in ms for each channel/trial
% stimOn is the relative stimulus onset time in ms for a given trial.
% blue and hemoAvg are the averages over all frames for each channel.

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

%construct path to data folder
opts.paradigm = 'SpatialDisc'; %this is the name of the paradigm. should be the same for all datasets.
opts.fPath = [cPath Animal filesep opts.paradigm filesep Rec]; %build file path
opts.frameRate = 30; %frame rate of individual channels. With alternating illumination this is the absolute frameRate / 2.
opts.preStim = 200 / opts.frameRate; %prestim duration in seconds.
opts.postStim = 200 / opts.frameRate; %poststim duration in seconds.
opts.animal = Animal; %animal ID
opts.rec = Rec; %name of current recording. usually the date of recording.
opts.fName = 'Frames'; %name of imaging data files.
opts.plotChans = false; %plot results for blue and violet channel when loading raw data. Will make the code break on the HPC.
opts.stimLine = 6; %analog line that contains stimulus trigger.
opts.barcodeLine = 7; %line for Bpod trialID.
opts.trigLine = [8 9]; %analog lines for blue and violet light triggers.
opts.memLimit = 5; % memory limit for video data in workspace in gigabyte. Use frame subsampling to stay within limit.
opts.loadRaw = true; %flag to use raw data instead of .mj2 files.
opts.verbosity = false; %flag to supress warnings from 'splitChannel / splitVideo' code.
opts.maxTrials = inf; %define a range of trials that should be used

%% check if .dat or .mj2 files are present. Use slower splitvideo code for latter case.
rawCheck = dir([opts.fPath filesep opts.fName '*dat']);
vidCheck = dir([opts.fPath filesep opts.fName '*mj2']);
analogCheck = dir([opts.fPath filesep 'Analog_*.dat']);

if size(rawCheck,1) == size(analogCheck,1)
elseif size(vidCheck,1) == size(analogCheck,1)
    opts.loadRaw = false;
    disp('Using .mj2 files instead of binary files.');
else
    error('Unequal number of imaging and analog data files. Aborted')
end

%% get some variables for analysis
disp(opts.fPath); disp(datestr(now));
load([opts.fPath filesep 'mask.mat'], 'mask'); %load mask for current session
load([opts.fPath filesep 'Vc.mat'], 'trials'); %load trials for current session
if length(trials) > opts.maxTrials
    fprintf('Using %d/%d recorded trials\n',opts.maxTrials, length(trials));
    trials = trials(1:opts.maxTrials);
end

frameCnt = repmat(round((opts.preStim + opts.postStim) * opts.frameRate), 1, length(trials));
opts.falseAlign = false(1,length(trials)); %this indicates that something unexpeced happened when loading raw data. Will always be true when asking for more frames as present in the raw files.

% get reference images for motion correction
[blueData,~,hemoData] = Widefield_SplitVideo(opts,trials(1));   
blueData = single(squeeze(blueData));
blueRef = fft2(median(blueData,3)); %blue reference for alignment
hemoData = single(squeeze(hemoData));
hemoRef = fft2(median(hemoData,3)); %violet reference for alignment

blueAvg = zeros([size(mask), length(trials)],'uint16'); %average for mean correction. Collect single session averages to verify correct channel separation.
hemoAvg = zeros([size(mask), length(trials)],'uint16');
    
blueData = arrayShrink(blueData(:,:,1),mask,'merge'); %get one compressed frame for size estimate
info = whos('blueData');
exptSize = info.bytes * sum(frameCnt) * 2 / 2^30; %expected size of complete data set in gb if all frames are present.
frameFrac = ceil(exptSize / opts.memLimit); %only use fraction of frames to keep memory usage under control
if frameFrac > 1 %adjust size estimate if averaging across frames
    mov = zeros(sum(~mask(:)),sum(floor(frameCnt/frameFrac)),2,'single'); %large mov matrix with reduced frame selection. Only use selected pixels later.
    expFrames = floor(round((opts.preStim + opts.postStim) * opts.frameRate)/frameFrac); %expected frames per trial after subsampling
    blueFrameTimes = NaN(1,sum(floor(frameCnt/frameFrac)));
    hemoFrameTimes = NaN(1,sum(floor(frameCnt/frameFrac)));
else
    mov = zeros(sum(~mask(:)),ceil(sum(frameCnt)),2,'single'); %large mov matrix. Only use selected pixels later.
    blueFrameTimes = NaN(1,ceil(sum(frameCnt/2)));
    hemoFrameTimes = NaN(1,ceil(sum(frameCnt/2)));
end

info = whos('mov');
redSize = info.bytes / 2^30;
fprintf(1,'MemLimit: %f gb ; Expected size: %f gb ; Selected fraction: %d ; Reduced size: %f gb \n', opts.memLimit,exptSize,frameFrac,redSize);
Cnt = 0;

%% perform image alignement for separate channels and collect data in mov matrix
for iTrials = 1:length(trials)
    
    [blueData,blueTimes,hemoData,hemoTimes,stimOn(iTrials),opts.falseAlign(iTrials)] = Widefield_SplitVideo(opts,trials(iTrials));
    
    if size(blueData,3) ~= size(hemoData,3)
        error(['Trial ' int2str(trials(iTrials)) ': Blue and hemo channels have uneven framecount'])
    end
    frameCnt(iTrials) = length(blueTimes); %adjust frameCnt, depending on how many frames were found
        
    if frameFrac > 1 %subselect frames if required
        [~,selInd] = sort(rand(frameCnt(iTrials),1)); % get random index for current trial
        if length(selInd) >= expFrames
            selInd = selInd(1:expFrames);
        end
        % only use subselection of data
        blueData = blueData(:,:,selInd); 
        hemoData = hemoData(:,:,selInd);
        blueTimes = blueTimes(selInd);
        hemoTimes = hemoTimes(selInd);
        
        %keep current frame selection for later
        frameIdx{iTrials} = selInd; 
    end
    
    %perform image alignment for both channels
    for iFrames = 1:size(blueData,3)
        [~, temp] = Widefield_dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10);
        blueData(:, :, iFrames) = abs(ifft2(temp));
        
        [~, temp] = Widefield_dftregistration(hemoRef, fft2(hemoData(:, :, iFrames)), 10);
        hemoData(:, :, iFrames) = abs(ifft2(temp));
    end
       
    % keep avg for each trial to check if channels were separated correctly
    blueAvg(:,:,iTrials) = median(blueData,3);
    hemoAvg(:,:,iTrials) = median(hemoData,3);
    
    %only use selected pixels from mask
    blueData = arrayShrink(blueData,mask,'merge'); 
    hemoData = arrayShrink(hemoData,mask,'merge');
    
    %keep timestamps for all frames
    blueFrameTimes(1,Cnt + (1:length(blueTimes))) = blueTimes;
    hemoFrameTimes(1,Cnt + (1:length(hemoTimes))) = hemoTimes;

    mov(:,Cnt + (1:size(blueData,2)),1) = blueData;    
    mov(:,Cnt + (1:size(blueData,2)),2) = hemoData;
    Cnt = Cnt + size(blueData,2);

    if rem(iTrials,100) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
end
clear blueData hemoData blueRef hemoRef

%save all averages to check for channel separation - this should be ok
% save([opts.fPath filesep 'allBlueAvg.mat'],'blueAvg');
% save([opts.fPath filesep 'allHemoAvg.mat'],'hemoAvg');
 
%compute average across all sessions for later data correction
blueAvg = mean(single(blueAvg),3);
hemoAvg = mean(single(hemoAvg),3);
% save([opts.fPath filesep 'blueAvg.mat'],'blueAvg');
% save([opts.fPath filesep 'hemoAvg.mat'],'hemoAvg');
% save([opts.fPath filesep 'frameIdx.mat'],'frameIdx'); %save index for sub selection

if size(mov,2) > Cnt(1) 
    fprintf(1, 'Counter = %d; mov size = (%d,%d)\n', Cnt(1),size(mov,1),size(mov,2));
    fprintf(1, 'Removing %d/%d entries from mov matrix\n', size(mov,2) - Cnt(1),size(mov,2));
    mov(:,Cnt(1)+1:end,:) = []; %remove unused entries from mov matrix
end
blueFrameTimes(isnan(blueFrameTimes)) = [];
hemoFrameTimes(isnan(hemoFrameTimes)) = [];

%% subtract and divide by basline - do this in steps to reduce memory load
blueAvg = arrayShrink(blueAvg,mask,'merge');
hemoAvg = arrayShrink(hemoAvg,mask,'merge');

ind = 0:200:size(mov,2);
for x = 1:length(ind)
    if x == length(ind)
        mov(:,ind(x)+1:end,1) = bsxfun(@minus, mov(:,ind(x)+1:end,1), blueAvg); % subtract blue mean 
        mov(:,ind(x)+1:end,1) = bsxfun(@rdivide, mov(:,ind(x)+1:end,1), blueAvg); % divide by blue mean 
    
        mov(:,ind(x)+1:end,2) = bsxfun(@minus, mov(:,ind(x)+1:end,2), hemoAvg); % subtract hemo mean 
        mov(:,ind(x)+1:end,2) = bsxfun(@rdivide, mov(:,ind(x)+1:end,2), hemoAvg); % divide by hemo mean 
    else
        mov(:,ind(x)+1:ind(x+1),1) = bsxfun(@minus, mov(:,ind(x)+1:ind(x+1),1), blueAvg); % subtract blue mean 
        mov(:,ind(x)+1:ind(x+1),1) = bsxfun(@rdivide, mov(:,ind(x)+1:ind(x+1),1), blueAvg); % divide blue mean 
       
        mov(:,ind(x)+1:ind(x+1),2) = bsxfun(@minus, mov(:,ind(x)+1:ind(x+1),2), hemoAvg); % subtract hemo mean
        mov(:,ind(x)+1:ind(x+1),2) = bsxfun(@rdivide, mov(:,ind(x)+1:ind(x+1),2), hemoAvg); % divide hemo mean 
    end
end
end