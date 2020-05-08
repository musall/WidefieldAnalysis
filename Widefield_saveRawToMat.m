function [frameCnt, blueFrameTimes, hemoFrameTimes, stimOn, blueAvg, hemoAvg, mask, opts] = Widefield_saveRawToMat(cPath, savePath, recType)
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

if ~exist('savePath', 'var') || isempty(savePath)
    savePath = cPath;
end

if ~exist('recType', 'var') || isempty(recType)
    recType = 'bpod';
end

%construct path to data folder
opts.fPath = cPath; %build file path
opts.fName = 'Frames'; %name of imaging data files.
opts.plotChans = false; %plot results for blue and violet channel when loading raw data. Will make the code break on the HPC.
opts.loadRaw = true; %flag to use raw data instead of .mj2 files.
opts.verbosity = false; %flag to supress warnings from 'splitChannel / splitVideo' code.
disp('===============');
disp(opts.fPath); disp(datestr(now));

if strcmpi(recType, 'bpod')
    disp('BpodImager settings');
    opts.stimLine = 6; %analog line that contains stimulus trigger.
    opts.trigLine = [8 9]; %analog lines for blue and violet light triggers.
elseif strcmpi(recType, 'widefield')
    disp('WidefieldImager settings');
    opts.stimLine = 3; %analog line that contains stimulus trigger.
    opts.trigLine = [6 7]; %analog lines for blue and violet light triggers.
end
disp('===============');

%% check if .dat or .mj2 files are present. Use slower splitvideo code for latter case.
rawCheck = dir([opts.fPath filesep opts.fName '*dat']);
vidCheck = dir([opts.fPath filesep opts.fName '*mj2']);
analogCheck = dir([opts.fPath filesep 'Analog_*.dat']);

if size(rawCheck,1) == size(analogCheck,1)
    fileCnt = size(rawCheck,1);
    for iFiles = 1 : fileCnt %go through files and check for trialNr. this is used for transformation to .mat files later.
        temp = textscan(rawCheck(iFiles).name,'%s%f%s','Delimiter','_');
        trials(iFiles) = temp{2};
    end
elseif size(vidCheck,1) == size(analogCheck,1)
    fileCnt = size(vidCheck,1);
    for iFiles = 1 : fileCnt %go through files and check for trialNr. this is used for transformation to .mat files later.
        temp = textscan(vidCheck(iFiles).name,'%s%f%s','Delimiter','_');
        trials(iFiles) = temp{2};
    end
    opts.loadRaw = false;
    disp('Using .mj2 files instead of binary files.');
else
    error('Unequal number of imaging and analog data files. Aborted')
end
trials = sort(trials);

%% get reference images for motion correction
if opts.loadRaw
    [blueData,~,hemoData] = Widefield_SplitChannels(opts,trials(1));
else
    [blueData,~,hemoData] = Widefield_SplitVideo(opts,trials(1));
end
blueData = single(squeeze(blueData));
blueRef = fft2(median(blueData,3)); %blue reference for alignment
save([savePath 'blueRef.mat'],'blueRef');
hemoData = single(squeeze(hemoData));
hemoRef = fft2(median(hemoData,3)); %violet reference for alignment
save([savePath 'hemoRef.mat'],'hemoRef');

blueAvg = zeros([size(blueData,1), size(blueData,2), fileCnt],'uint16'); %average for mean correction. Collect single session averages to verify correct channel separation.
hemoAvg = zeros([size(blueData,1), size(blueData,2), fileCnt],'uint16');
blueFrameTimes = cell(1, fileCnt);
hemoFrameTimes = cell(1, fileCnt);

%% perform image alignement for separate channels and collect data in mov matrix
for iTrials = 1:fileCnt
    
    if opts.loadRaw
        [blueData,blueTimes,hemoData,hemoTimes] = Widefield_SplitChannels(opts,trials(iTrials));
    else
        [blueData,blueTimes,hemoData,hemoTimes] = Widefield_SplitVideo(opts,trials(iTrials));
    end

    if size(blueData,3) ~= size(hemoData,3)
        error(['Trial ' int2str(trials(iTrials)) ': Blue and hemo channels have uneven framecount'])
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
    hemoAvg(:,:,iTrials) = median(blueData,3);
    
    %keep timestamps for all frames
    blueFrameTimes{iTrials} = blueTimes;
    hemoFrameTimes{iTrials} = hemoTimes;

    if rem(iTrials,25) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
    
    save([savePath filesep strcat('blue_trial_',num2str(trials(iTrials)),'.mat')],'blueData', 'blueTimes');
    save([savePath filesep strcat('hemo_trial_',num2str(trials(iTrials)),'.mat')],'hemoData', 'hemoTimes');
end
clear blueData hemoData blueTimes hemoTiems blueRef hemoRef 

save([savePath 'trials.mat'],'trials'); %save trials so order of analysis is consistent

%save frametimes for blue/hemo trials
save([savePath 'blueFrameTimes.mat'],'blueFrameTimes');
save([savePath 'hemoFrameTimes.mat'],'hemoFrameTimes');

%save averages in case you need them later
save([savePath 'blueAvg.mat'],'blueAvg');
save([savePath 'hemoAvg.mat'],'hemoAvg');

%take average over all trials for subsequent mean correction
blueAvg = mean(single(blueAvg),3);
hemoAvg = mean(single(hemoAvg),3);
            
%% subtract and divide by baseline - do this per trial to reduce memory load
for iTrials = 1:length(trials)
    load([savePath filesep strcat('blue_trial_',num2str(trials(iTrials)),'.mat')])
    load([savePath filesep strcat('hemo_trial_',num2str(trials(iTrials)),'.mat')])

    blueDataMat=(single(blueData)-blueAvg)./blueAvg;
    save([savePath filesep strcat('blue_trial_',num2str(trials(iTrials)),'.mat')],'blueDataMat');
    hemoDataMat=(single(hemoData)-hemoAvg)./hemoAvg;
    save([savePath filesep strcat('hemo_trial_',num2str(trials(iTrials)),'.mat')],'hemoDataMat');
end
