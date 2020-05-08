function [frameCnt, blueFrameTimes, hemoFrameTimes, stimOn, blueAvg, hemoAvg, mask, opts] = Widefield_saveToMat(cPath, savepath)
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

if ~exist(savePath) || isempty(savePath)
    savePath = cPath;
end

%construct path to data folder
opts.fPath = cPath; %build file path
opts.frameRate = 30; %frame rate of individual channels. With alternating illumination this is the absolute frameRate / 2.
opts.preStim = 200 / opts.frameRate; %prestim duration in seconds.
opts.postStim = 200 / opts.frameRate; %poststim duration in seconds.
opts.fName = 'Frames'; %name of imaging data files.
opts.plotChans = false; %plot results for blue and violet channel when loading raw data. Will make the code break on the HPC.
opts.stimLine = 6; %analog line that contains stimulus trigger.
opts.barcodeLine = 7; %line for Bpod trialID.
opts.trigLine = [8 9]; %analog lines for blue and violet light triggers.
opts.memLimit = 210; % memory limit for video data in workspace in gigabyte. Use frame subsampling to stay within limit.
opts.loadRaw = true; %flag to use raw data instead of .mj2 files.
opts.verbosity = false; %flag to supress warnings from 'splitChannel / splitVideo' code.

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


frameCnt = repmat(round((opts.preStim + opts.postStim) * opts.frameRate), 1, length(trials));
opts.falseAlign = false(1,length(trials)); %this indicates that something unexpeced happened when loading raw data. Will always be true when asking for more frames as present in the raw files.

% get reference images for motion correction
[blueData,~,hemoData] = Widefield_SplitVideo(opts,trials(1));   
blueData = single(squeeze(blueData));
blueRef = fft2(median(blueData,3)); %blue reference for alignment
hemoData = single(squeeze(hemoData));
hemoRef = fft2(median(hemoData,3)); %violet reference for alignment
    
blueData = arrayShrink(blueData(:,:,1),mask,'merge'); %get one compressed frame for size estimate
info = whos('blueData');
exptSize = info.bytes * sum(frameCnt) * 2 / 2^30; %expected size of complete data set in gb if all frames are present.
frameFrac = ceil(exptSize / opts.memLimit); %only use fraction of frames to keep memory usage under control
if frameFrac > 1 %adjust size estimate if averaging across frames
    blueFrameTimes = NaN(1,sum(floor(frameCnt/frameFrac)));
    hemoFrameTimes = NaN(1,sum(floor(frameCnt/frameFrac)));
else
    blueFrameTimes = NaN(1,ceil(sum(frameCnt/2)));
    hemoFrameTimes = NaN(1,ceil(sum(frameCnt/2)));
end

Cnt = 0;

blueAvg = zeros([size(mask), length(trials)],'uint16'); %average for mean correction. Collect single session averages to verify correct channel separation.
hemoAvg = zeros([size(mask), length(trials)],'uint16');

%% perform image alignement for separate channels and collect data in mov matrix
for iTrials = 1:length(trials)
    if exist([savepath filesep strcat('blue_trial_',num2str(trials(iTrials)),'.mat')],'file'), continue; end
    try
        [blueData,blueTimes,hemoData,hemoTimes,stimOn(iTrials),opts.falseAlign(iTrials)] = Widefield_SplitVideo(opts,trials(iTrials));
    catch
        continue;
    end
    if size(blueData,3) ~= size(hemoData,3)
        error(['Trial ' int2str(trials(iTrials)) ': Blue and hemo channels have uneven framecount'])
    end
    frameCnt(iTrials) = length(blueTimes); %adjust frameCnt, depending on how many frames were found
            
    %perform image alignment for both channels
    for iFrames = 1:size(blueData,3)
        [~, temp] = Widefield_dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10);
        blueData(:, :, iFrames) = abs(ifft2(temp));
        
        [~, temp] = Widefield_dftregistration(hemoRef, fft2(hemoData(:, :, iFrames)), 10);
        hemoData(:, :, iFrames) = abs(ifft2(temp));
    end
       
    % keep avg for each trial to check if channels were separated correctly
    blueAvg(:,:,iTrials) = median(blueData,3);
    hemoAvg(:,:,iTrials) = median(blueData,3);
    
    
    %keep timestamps for all frames
    blueFrameTimes(1,Cnt + (1:length(blueTimes))) = blueTimes;
    hemoFrameTimes(1,Cnt + (1:length(hemoTimes))) = hemoTimes;

    Cnt = Cnt + size(blueData,3);

    if rem(iTrials,100) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
    
    save([savepath filesep strcat('blue_trial_',num2str(trials(iTrials)),'.mat')],'blueData');
    save([savepath filesep strcat('hemo_trial_',num2str(trials(iTrials)),'.mat')],'hemoData');
end
clear blueData hemoData blueRef hemoRef 

blueFrameTimes(isnan(blueFrameTimes)) = [];
hemoFrameTimes(isnan(hemoFrameTimes)) = [];
save([savepath filesep 'Avg_mats.mat'],'blueAvg','hemoAvg');

blueAvg = mean(single(blueAvg),3);
hemoAvg = mean(single(hemoAvg),3);

% %% subtract and divide by baseline - do this per trial to reduce memory load
for iTrials = 1:length(trials)
    try
        load([savepath filesep strcat('blue_trial_',num2str(trials(iTrials)),'.mat')])
        load([savepath filesep strcat('hemo_trial_',num2str(trials(iTrials)),'.mat')])
    catch 
        continue
    end
    blueDataMat=(blueData-blueAvg)./blueAvg;
    save([savepath filesep strcat('blue_trial_',num2str(trials(iTrials)),'.mat')],'blueDataMat');
    hemoDataMat=(hemoData-hemoAvg)./hemoAvg;
    save([savepath filesep strcat('hemo_trial_',num2str(trials(iTrials)),'.mat')],'hemoDataMat');
end
end
