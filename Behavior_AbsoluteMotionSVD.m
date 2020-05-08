function Behavior_AbsoluteMotionSVD(fPath, Animal, Rec, bhvOpts, movType)
% code to convert behavioral videos from timestamp (.mat) and compressed
% movie (.mj2) files. This codes computes an SVD representation of the
% absolute first derivative of video data to compute change in animal motion.


nrCams = 2;
if ~exist('bhvOpts','var') || isempty(bhvOpts)
    bhvOpts.nSVD = 500;
    bhvOpts.maxFrameCnt = 500; %max number of frames per trial. Earlier frames will be removed from analysis.
    bhvOpts.memLimit = 32; % memory limit for video data in workspace in gigabyte. Determines how many segments can be computed at once.
    bhvOpts.segSize = 4; %segment size relative to full video. Minimum is 2. Will create segSize^2 segments and create svd for each. Second SVD combines all segments into one V/U.
    bhvOpts.targRate = 30; %rate of widefield imaging in Hz. For combined SVD, video data is resampled to match targRate.
    bhvOpts.preStimDur = 3; %baseline duration for trial-aligned video data
    bhvOpts.postStimDur = 115/bhvOpts.targRate; % Duration of trial after stimulus onset in seconds
    bhvOpts.framesPerTrial = round((bhvOpts.preStimDur + bhvOpts.postStimDur) * bhvOpts.targRate); %amount of frames per trial in trial-aligned video data
    bhvOpts.trialDur = bhvOpts.framesPerTrial / bhvOpts.targRate; %amount of frames per trial in trial-aligned video data
    bhvOpts.smth = 5; % filter length when smoothing video data before combined svd
end

if isfield(bhvOpts,'fPath') % if bhvOpts has fPath already, use it from there. Otherwise use current input.
    fPath = bhvOpts.fPath;
    if fPath(end) ~= filesep
        fPath = [fPath filesep];
    end
else
    if fPath(end) ~= filesep
        fPath = [fPath filesep];
    end
    fPath = [fPath Animal filesep 'SpatialDisc' filesep Rec filesep 'BehaviorVideo' filesep];
end
if ~exist('movType','var') || isempty(movType)
    movType = '.mj2';
end

%% load behavior data
cPath = [fileparts(fPath(1:end-1)) filesep];
[~, cFile] = fileparts(fPath(1:end-1)); %behavioral file name must match name of the data folder.
load([cPath cFile]);

%% loop through available webcams
for iCams = 1:nrCams
    
    timeFiles = dir([fPath '*frameTimes*_' int2str(iCams) '.mat']); %all files for current cam based on integer before .mat
    movieFiles = dir([fPath '*Video*_' int2str(iCams) movType]); %all files for current cam based on integer before .mat
    trials = 1:length(movieFiles);
    
    if length(movieFiles) > 1 %don't do conversion if less than 25 movie files are found - probably something was not deleted correctly
        % check frameTimes to asses total number of frames in experiment
        frameCnt = 0;
        timeStamps = cell(1,length({timeFiles.name}));
        for iFiles = 1:length(timeFiles)
            cFile = [fPath timeFiles(iFiles).name];
            load(cFile)
            if size(frameTimes,1) > bhvOpts.maxFrameCnt
                timeStamps{iFiles} = frameTimes(end-bhvOpts.maxFrameCnt+1:end);
            else
                timeStamps{iFiles} = frameTimes;
            end
            frameCnt(iFiles) = size(timeStamps{iFiles},1);  %nr of frames per trial
        end
        totalFrameTimes = cat(1,timeStamps{:});
        save([fPath 'motionSVD_Cam' int2str(iCams) '-frameTimes.mat'],'totalFrameTimes','frameCnt');
        
        %get single frame to assess video size
        cFile = [fPath movieFiles(1).name];
        v = VideoReader(cFile);
        singleFrame = single(arrayResize(readFrame(v),2)); clear v
        if size(singleFrame,3) > 1; singleFrame = singleFrame(:,:,1); end
        
        blockSize = ceil(size(singleFrame)/bhvOpts.segSize);
        ind = im2col(reshape(1:numel(singleFrame),size(singleFrame)),blockSize,'distinct'); %build index for diffent image blocks
        save([fPath 'motion_segInd' int2str(iCams) '.mat'],'ind'); %save index for later reconstruction
        
        %% determine how often to load raw data given available memory
        info = whos('singleFrame');
        exptSize = (info.bytes * sum(frameCnt) / 2^30) / size(ind,2); %expected size in gb for each video segment
        segPerRun = floor(bhvOpts.memLimit / exptSize); %number of segments per run
        segRuns = ceil(size(ind,2) / segPerRun); %required amount of loading raw data to stay within memory constraints
        fprintf(1, 'Cam%d: Expected datasize: %fgb. Reload data %d times.\n',iCams, exptSize, segRuns);
        segCnt = 0;
        
        %% compute svd for each video segment
        for iSegs = 1 : segRuns
            if iSegs == segRuns %for last run, check if segPerRun should be reduced
                segPerRun = size(ind,2) - segCnt;
            end
            for x = 1 : segPerRun
                mov{x} = zeros(sum(frameCnt),sum(ind(:,segCnt + x) <= numel(singleFrame)),'single'); %pre-allocate mov array to compute svd later. Allows for different segment sizes if index is larger as a single frame.
            end
            
            Cnt = 0;
            for iFiles = 1:length(movieFiles)
                cFile = [fPath movieFiles(iFiles).name];
                rawData = squeeze(importdata(cFile));
                if size(rawData,3) == 3; rawData = squeeze(rawData(:,:,1,:));end
                
                if size(rawData,3) > bhvOpts.maxFrameCnt
                    rawData = rawData(:,:,end-bhvOpts.maxFrameCnt+1:end);
                end
                rawData = cat(3,rawData(:,:,1),rawData); %duplicate first frame
                rawData = abs(diff(rawData,[],3));
                
                %% merge data in mov array for svd compression
                rawData = arrayResize(rawData,2);
                rawData = reshape(rawData, [], size(rawData,3)); %reshape data
                
                for x = 1 : segPerRun
                    cInd = ind(:,segCnt + x);
                    cInd(cInd > size(rawData,1)) = []; %ensure index does not go beyond available pixels
                    mov{x}(Cnt + (1:size(rawData,2)), :) =  rawData(cInd,:)';  %get current segment
                end
                Cnt = Cnt + size(rawData,2);
                
                if rem(iFiles,50) == 0
                    fprintf(1, 'Cam%d: Loaded file %d/%d\n',iCams, iFiles,length(trials));
                end
            end
            
            if Cnt ~= size(mov{1},1)
                error('Something is wrong with mov array. Cnt ~= size(mov,1).')
            end
            
            %% compute svd for current segments
            for x = 1 : segPerRun
                tic; fprintf(1, 'Computing SVD - Cam %d, Segment %d\n',iCams,segCnt + x);
                [V,S,U]  = fsvd(mov{x}, bhvOpts.nSVD);
                V = V * S;
                V = V(:,1:bhvOpts.nSVD); %get request nr of dimensions
                U = U(:,1:bhvOpts.nSVD)'; %the resulting U is nSVD x regs
                save([fPath 'motionSVD_Cam' int2str(iCams) '-Seg' int2str(segCnt + x) '.mat'],'V','U');
                toc
            end
            clear mov V U frameTimes
            segCnt = segCnt + segPerRun;
            
        end
    end
end

%% load svd data and merge different segments together
allV = cell(1,nrCams);
allTimes = cell(1,nrCams);
for iCams = 1:nrCams %get segments and frameTimes for each cam
    load([fPath 'motionSVD_Cam' int2str(iCams) '-frameTimes.mat'],'totalFrameTimes','frameCnt');
    load([fPath 'motion_segInd' int2str(iCams) '.mat'],'ind');
    
    % check if timestamps are shifted by an hour. Apparently that can happen.
    timeCheck = (SessionData.TrialStartTime(1)*86400) - (totalFrameTimes(1) * 86400); %time difference between first acquired frame and onset of first trial
    if timeCheck > 3540 && timeCheck < 3660 %timeshift by one hour (+- 60seconds)
        warning('Behavioral and video timestamps are shifted by 1h. Will adjust timestamps in video data accordingly.')
        totalFrameTimes = totalFrameTimes + (1/24); %add one hour
    elseif timeCheck > 30
        error('Something wrong with timestamps in behavior and video data. Time difference is larger as 30 seconds.')
    end
    
    allTimes{iCams} = totalFrameTimes; %get timestamps for current cam
    allFrameCnt{iCams} = frameCnt;
    bhvOpts.frameRate(iCams) = round(1/median(diff(allTimes{iCams} * 86400))); %determine frameRate for current cam
    
    % load compressed video segments
    allV{iCams} = cell(1,size(ind,2));
    for iSegs = 1 : size(ind,2)
        data = load([fPath 'motionSVD_Cam' int2str(iCams) '-Seg' int2str(iSegs) '.mat'],'V');
        allV{iCams}{iSegs} = data.V; clear data
    end
    allV{iCams} = cat(2,allV{iCams}{:}); % combine in on larger array
end

%% go through each trial and build combined dataset for each cam
for x = 1 : nrCams
    alignV{x} = cell(1,size(SessionData.RawEvents.Trial,2));
end

for iTrials = 1:size(SessionData.RawEvents.Trial,2)
    try
        stimTime(iTrials) = SessionData.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
    catch
        stimTime(iTrials) = NaN;
    end
    
    for iCams = 1:nrCams
        try
            trialOn = (SessionData.TrialStartTime(iTrials) * 86400) + (stimTime(iTrials) - bhvOpts.preStimDur); % trial onset times
            cTimes = allTimes{iCams} * 86400 - trialOn; % get frame times for current cam
            trialOn = find(cTimes > 0,1); %trial on in frames
            trialOff = find(cTimes > bhvOpts.trialDur,1); %trial off in frames
            
            trialTime = cTimes(trialOn:trialOff) - cTimes(trialOn);
            idx = trialTime < bhvOpts.trialDur; %don't use late frames
            trialTime = trialTime(idx);
            timeLeft = bhvOpts.trialDur - trialTime(end); %check if there is missing time at the end of a trial
            
            cData = allV{iCams}(trialOn:trialOff,:); %get data based on time vector matching duration of a trial
            cData = cData(idx, :); %don't use late frames
            
            if (timeLeft < bhvOpts.trialDur * 0.9) && (timeLeft > 0) %if there is some time missing to make a whole trial
                addTime = trialTime(end) + (1/bhvOpts.targRate : 1/bhvOpts.targRate : timeLeft + 1/bhvOpts.targRate); %add some dummy times to make complete trial
                trialTime = [trialTime' addTime];
                cData = [cData; repmat(cData(1,:),length(addTime),1)];
            end
            
            if any(diff(trialTime) > 0.5) %dont use current trial if frames are no more than 0.5s apart
                error;
            end
            
            %create time vector to capture potential irregular sampling
            timePad = 1/bhvOpts.targRate:1/bhvOpts.targRate:1/bhvOpts.targRate*21;
            trialTime = trialTime'+timePad(end)+1/bhvOpts.frameRate(iCams);
            trialTime = [timePad trialTime' timePad+trialTime(end)];
            cData = [repmat(cData(1,:),21,1); cData; repmat(cData(end,:),21,1)]; %add some padding on both sides to avoid edge effects when resampling
            cData = resample(double(cData), trialTime, bhvOpts.targRate); %resample to set target framerate
            cData = conv2(cData',ones(1,bhvOpts.smth)/bhvOpts.smth,'same')'; %smooth trace with moving average of 'smth' points
            alignV{iCams}{iTrials} = cData(22 : 21 + bhvOpts.framesPerTrial,:); %remove padds and select trial-relevant frames
        catch
            alignV{iCams}{iTrials} = NaN(bhvOpts.framesPerTrial,size(allV{iCams},2),'single'); %dont use this trial
        end
    end
end
clear allV

% combined V segments from all cams
for x = 1 : nrCams
    alignV{x} = cat(1,alignV{x}{:}); %combined V matrix
end
alignV = cat(2,alignV{:}); %combined V matrix
camIdx = ~isnan(mean(alignV,2)); %index for trials where bhv frames were available
meanTrace = nanmean(alignV(camIdx,:),2); %get mean change in both cameras

% create second SVD
tic; disp('Compute combined SVD.');
[V,S,vidU]  = fsvd(alignV(camIdx,:), bhvOpts.nSVD); clear alignV
V = V * S;
V = V(:,1:bhvOpts.nSVD); %get request nr of dimensions
vidU = vidU(:,1:bhvOpts.nSVD)'; %the resulting U is nSVD x regs

%zscore video regressor, save mean and std for later reconstruction
meanV = mean(V);
stdV = std(V);
V = bsxfun(@minus,V,meanV);
V = bsxfun(@rdivide,V,stdV);

% rebuild
vidV = NaN(length(camIdx),bhvOpts.nSVD,'single');
vidV(camIdx,:) = V; clear V
save([fPath 'motionSVD_CombinedSegments.mat'], 'vidV', 'vidU', 'meanV' ,'stdV','allTimes');
save([fPath 'motionSVD_meanTrace.mat'], 'meanTrace');
save([fPath 'motionbhvOpts.mat'], 'bhvOpts');
disp('Done.'); toc;

end