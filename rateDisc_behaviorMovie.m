% rateDisc_behaviorMovie

load allenDorsalMapSM
allenMask = dorsalMaps.allenMask;
sRate = 30;           % Sampling rate of imaging in Hz
preStimDur = floor(2 * sRate) / sRate; % Duration of trial before lever grab in seconds
postStimDur = floor(3 *sRate) / sRate; % Duration of trial after lever grab onset in seconds
frames = round((preStimDur + postStimDur) * sRate); %nr of frames per trial
trialDur = (frames * (1/sRate)); %duration of trial in seconds

%% get imaging data for example movie
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';
animal = 'mSM66';
cRec = '14-Jun-2018';
% cPath = 'X:\smusall\BpodImager\Animals\';
% animal = 'mSM78';
% cRec = '06-Sep-2018';

fPath = [cPath animal filesep 'SpatialDisc' filesep cRec filesep];
load([fPath 'interpVc'],'Vc','frames');
load([fPath 'regData'],'trialIdx');
load([fPath 'Vc'],'U','bTrials');
load([fPath 'mask'],'mask');
load([fPath 'opts2'],'opts');
trialIdx = unique(ceil(find(~trialIdx)/frames)); %rebuild as single trial index (instead of frames)
bTrials = bTrials(trialIdx);

%% align data to allen
U = alignAllenTransIm(U, opts.transParams);
U = arrayShrink(U(:,1:size(allenMask,2),:),allenMask,'merge');
trialAvg = nanmean(U,1) * Vc;
trialAvg = reshape(trialAvg,frames,[]);
[a,b] = sort(nanmean(trialAvg,1),'ascend'); %this is optional to isolate high-activity frames. Not really important though.
% trialList = [100 b(end)];
trialList = 1:5; %trials to use

%% make frame sequence
Vc = reshape(Vc,size(Vc,1),frames,[]);
temp = smoothCol(reshape(Vc(:,:,trialList),size(Vc,1),[]),2,'gauss',[],2);
cTrial = U * temp;
cTrial = arrayShrink(cTrial,allenMask,'split');

%% get timestamps and behavioral data
bhvFile = dir([fPath filesep animal '_*.mat']);
load([fPath bhvFile(1).name],'SessionData'); %load behavior data
SessionData.TrialStartTime = SessionData.TrialStartTime * 86400; %convert trailstart timestamps to seconds
bhv = selectBehaviorTrials(SessionData,bTrials); %only use completed trials that are in the Vc dataset

load([fPath 'BehaviorVideo' filesep 'FilteredPupil.mat'], 'pTime', 'fPupil', 'sPupil', 'whisker', 'faceM', 'bodyM', 'nose', 'bTime'); %load pupil data
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

%% get behavioral movies and align for each trial - this comes from the rateDisc_RegressModel code
% mPath = 'X:\smusall\BehaviorVideo\Animals\'; %path to behavioral videos
% fmPath = [mPath animal filesep cRec filesep];
fmPath = [fPath 'BehaviorVideo' filesep];

Cnt = 0;
clear imaging video1 video2 eyeTrace
for iTrials = trialList
    Cnt = Cnt + 1;
    
    % imaging data
    imaging{Cnt} = cTrial(:,:,(Cnt-1)*frames + (1 : frames)); % bring to same format as behavioral video
    imaging{Cnt} = cat(3,imaging{Cnt},NaN(size(cTrial,1),size(cTrial,2),2));
    
    % align to imaging data from model
    leverTimes = [reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(bhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
    
    stimGrab = leverTimes(find(leverTimes == bhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
    trialOn = bhv.TrialStartTime(iTrials) + (stimGrab - preStimDur);
    
    % cam 1
    cFile = dir([fmPath animal '*_Video_' num2str(bTrials(iTrials), '%04i') '_1.mp4']);
    rawData = squeeze(importdata([fmPath cFile.name]));
    rawData = rot90(squeeze(rawData(:,:,1,:)),3);
    if size(rawData,3) > 500
        rawData = rawData(:,:,end-500+1:end);
    end
    trialTime = pTime{bTrials(iTrials)} - trialOn;
    rejIdx = find(trialTime < trialDur); %don't use late frames
    video1{Cnt} = rawData(:,:,rejIdx(end-frames+1:end));
    video1{Cnt} = cat(3,video1{Cnt},NaN(size(rawData,1),size(rawData,2),2));
    
    % cam 2
    cFile = dir([fmPath animal '*_Video_' num2str(bTrials(iTrials), '%04i') '_2.mp4']);
    rawData = squeeze(importdata([fmPath cFile.name]));
    rawData = rot90(squeeze(rawData(:,:,1,:)),3);
    if size(rawData,3) > 500
        rawData = rawData(:,:,end-500+1:end);
    end
    trialTime = bTime{bTrials(iTrials)} - trialOn;
    rejIdx = find(trialTime < trialDur); %don't use late frames
    video2{Cnt} = rawData(:,:,rejIdx(end-frames+1:end));
    video2{Cnt} = cat(3,video2{Cnt},NaN(size(rawData,1),size(rawData,2),2));
    
    % eyeTrace
    cFile = dir([fmPath 'eyeTrace_' num2str(bTrials(iTrials)) '.mj2']);
    rawData = squeeze(importdata([fmPath cFile.name]));
    rawData = rot90(rawData,3);
    if size(rawData,3) > 500
        rawData = rawData(:,:,end-500+1:end);
    end
    trialTime = pTime{bTrials(iTrials)} - trialOn;
    rejIdx = find(trialTime < trialDur); %don't use late frames
    eyeTrace{Cnt} = rawData(:,:,rejIdx(end-frames+1:end));
    eyeTrace{Cnt} = cat(3,eyeTrace{Cnt},NaN(size(rawData,1),size(rawData,2),2));
    
end
imaging = cat(3,imaging{:});
video1 = cat(3,video1{:});
video2 = cat(3,video2{:});
eyeTrace = cat(3,eyeTrace{:});

%% make single trial movies
wfRange = [-0.1 0.05];
Widefield_SaveToAvi(imaging, 'imgExample', 15, 'colormap_blueblackred', wfRange, fPath, 4)

bhvRange = [0 255];
Widefield_SaveToAvi(video1, 'faceExample', 15, 'gray', bhvRange, fPath, 4)
Widefield_SaveToAvi(video2, 'bodyExample', 15, 'gray', bhvRange, fPath, 4)
Widefield_SaveToAvi(eyeTrace, 'eyeExample', 15, 'gray', bhvRange, fPath, 4)

%% run LEAP - facecam
vidSize = size(video1);
cVideo = reshape(video1, vidSize(1), vidSize(2), 1, vidSize(end));
cVideo = imresize(cVideo,0.5);
cVideo = cVideo(:,:,1,~(cVideo(1,1,1,:)==0));
if exist([fPath 'faceCam_exampleVideo.h5\'], 'file')
    delete([fPath 'faceCam_exampleVideo.h5\']);
end
    
h5create([fPath 'faceCam_exampleVideo.h5\'],'/box', size(cVideo),'ChunkSize', size(cVideo),'Datatype','uint8')
h5write([fPath 'faceCam_exampleVideo.h5\'], '/box', uint8(cVideo));
%label_joints %call LEAP if needed
modelPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\mSM66\SpatialDisc\14-Jun-2018\models\190616_170820-n=462\final_model.h5';
preds = predict_box(cVideo, modelPath);
preds = preds.positions_pred;
skeleton = load('\\grid-hs\churchland_nlsas_data\data\BpodImager\Skeletons\rotateFaceSkeleton.mat');

h = figure;
v = VideoWriter([fPath 'faceLabelsExampleSmall.mp4'], 'MPEG-4'); %save as compressed video file
v.Quality = 100;
v.FrameRate = 15;
open(v);
frameRange = [1:150 451:750];

% for iFrames = 1 : size(cVideo,4)
for iFrames = frameRange
    
    imshow(cVideo(:,:,1,iFrames),[0 255]); hold on;
    plot(preds(1:6,1,iFrames),preds(1:6,2,iFrames),'b.','MarkerSize',30)
    plot(preds(7:10,1,iFrames),preds(7:10,2,iFrames),'r.','MarkerSize',30)
    plot(preds(12:13,1,iFrames),preds(12:13,2,iFrames),'k.','MarkerSize',30)
    writeVideo(v,getframe(gca)); hold off;

    
    if rem(iFrames,frames) == 0
        imshow(cVideo(:,:,1,iFrames),[255 256]);
        for x = 1 : 2
            writeVideo(v,getframe(gca)); hold off;
        end
    end
end
close(v);


%% run LEAP - bodycam
vidSize = size(video2);
cVideo = reshape(video2, vidSize(1), vidSize(2), 1, vidSize(end));
cVideo = imresize(cVideo,0.5);
cVideo = cVideo(:,:,1,~(cVideo(1,1,1,:)==0));
cVideo = cVideo(:,:,1,frameRange);
if exist([fPath 'bodyCam_exampleVideoSmall.h5\'], 'file')
    delete([fPath 'bodyCam_exampleVideoSmall.h5\']);
end
h5create([fPath 'bodyCam_exampleVideoSmall.h5\'],'/box', size(cVideo),'ChunkSize', size(cVideo),'Datatype','uint8')
h5write([fPath 'bodyCam_exampleVideoSmall.h5\'], '/box', uint8(cVideo));

%label_joints %call LEAP if needed
modelPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\mSM66\SpatialDisc\14-Jun-2018\models\190618_003410-n=108\final_model.h5';
preds = predict_box(cVideo, modelPath);
preds = preds.positions_pred;
skeleton = load('\\grid-hs\churchland_nlsas_data\data\BpodImager\Skeletons\rotateBodySkeleton.mat');

h = figure;
v = VideoWriter([fPath 'bodyLabelsExampleSmall.mp4'], 'MPEG-4'); %save as compressed video file
v.Quality = 100;
v.FrameRate = 15;
open(v);

for iFrames = 1 : size(cVideo,4)
    
    imshow(cVideo(:,:,1,iFrames),[0 255]); hold on;
    plot(preds(1:2,1,iFrames),preds(1:2,2,iFrames),'k.','MarkerSize',30)
    plot(preds(3:5,1,iFrames),preds(3:5,2,iFrames),'r.','MarkerSize',30)
    plot(preds(6,1,iFrames),preds(6,2,iFrames),'g.','MarkerSize',30)
    plot(preds(8,1,iFrames),preds(8,2,iFrames),'b.','MarkerSize',30)
    writeVideo(v,getframe(gca)); hold off;

    
    if rem(iFrames,frames) == 0
        imshow(cVideo(:,:,1,iFrames),[255 256]);
        for x = 1 : 2
            writeVideo(v,getframe(gca)); hold off;
        end
    end
end
close(v);


%% apply tracking for behavior movement video
eyeVars.orient = NaN(size(eyeTrace,3),1);
eyeVars.center = NaN(size(eyeTrace,3),2);
eyeVars.axes = NaN(size(eyeTrace,3),2);
eyeVars.solidity = NaN(size(eyeTrace,3),1);
highTresh = 0.25;
lowThresh = 0;

% get first estimate
eyeTrace = mat2gray(eyeTrace);
for iFrames = 1:size(eyeTrace,3)
    if eyeTrace(1,1,iFrames) ~= 0
        temp = Behavior_EyeCheck(eyeTrace(:,:,iFrames),highTresh,lowThresh,false); %get pupil data
        eyeVars.axes(iFrames,:) = temp.axes; %circle axis
        eyeVars.center(iFrames,:) = temp.center; %circle center
        eyeVars.orient(iFrames) = temp.orient; %circle center
    else
        eyeVars.axes(iFrames,:) = NaN(1,2); %circle axis
        eyeVars.center(iFrames,:) = NaN(1,2); %circle center
        eyeVars.orient(iFrames) = NaN; %circle center
    end
end

% find bad pupil estimates and try to repair.
trace = mean(eyeVars.axes,2);
if sum(trace == 0) ~= length(trace)
    trace(trace == 0) = NaN;
    trace = fillgaps(trace, 100);
    smoothTrace = smooth(trace,10,'rlowess'); %use smooth pupil trace as reference for outliers
else
    smoothTrace = trace;
end
smoothDiff = abs(mean(eyeVars.axes,2)-smoothTrace); %difference from smoothed pupil trace
smoothStd = nanstd(smoothDiff);
smoothDiff = smoothDiff ./ smoothStd; %put to STDs
% idx = find([0;diff(zscore(mean(eyeVars.axes,2)))]<-.5 | [0;diff(zscore(mean(eyeVars.axes,2)))]>.5 | smoothDiff > 0.5); %potentially bad results. Try to do better by finding result that is closest to average.
idx = 1:size(eyeVars.orient,1);
idx(isnan(eyeVars.orient)) = [];

for iFrames = 1:length(idx)
    clear temp tDiff
    Cnt = 0;
    for x = [-6 -2 2 6] %try eroding or fusing patches to get better pupil estimate
        Cnt = Cnt + 1;
        temp{Cnt} = Behavior_EyeCheck(eyeTrace(:,:,idx(iFrames)),highTresh,lowThresh,false,true,x); %get pupil data. Force individual areas to fuse together.
        tDiff(Cnt) = abs(smoothTrace(idx(iFrames)) - mean(temp{Cnt}.axes,2)) ./ smoothStd; %find modification that gets estimate closer to smoothed trace
    end
    [~,minDiff] = min(tDiff);
    if tDiff(minDiff) < abs(trace(idx(iFrames))-smoothTrace(idx(iFrames)))
        eyeVars.center(idx(iFrames),:) = temp{minDiff}.center; %circle center
        eyeVars.axes(idx(iFrames),:) = temp{minDiff}.axes; %circle axis
        eyeVars.orient(idx(iFrames)) = temp{minDiff}.orient; %circle orientation
        eyeVars.solidity(idx(iFrames)) = temp{minDiff}.solidity; %circle solidity
    end
end

% make pupil ellipse movie
imgScale = 5;
phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);
R = [ cos(theta)   sin(theta)
    -sin(theta)   cos(theta)];

figure
for iFrames = 1:size(eyeVars.orient,1)
    theta = pi*eyeVars.orient(iFrames)/180;
    xy = [(eyeVars.axes(iFrames,1)/2)*cosphi; (eyeVars.axes(iFrames,2)/2)*sinphi];
    xy = R*xy;
    
    x = xy(1,:) + eyeVars.center(iFrames,1);
    y = xy(2,:) + eyeVars.center(iFrames,2);
    pic = eyeTrace(:,:,iFrames);
    pic = imresize(pic, imgScale);
    
    if ~any(eyeVars.axes(iFrames,:) < 8)
        
        a = interp(x*imgScale,imgScale*2);
        b = interp(y*imgScale,imgScale*2);
        a = round(a); b = round(b);
        
        idx = a < size(pic,1) | b < size(pic,2);
        a = a(idx);
        b = b(idx);
        idx = sub2ind(size(pic), b, a);
        pic(idx) = 1;
    end
    
    %     imagesc(pic); colormap gray; axis square
    %     pause
    newEyeTrace(:,:,iFrames) = pic;
end

%% make single trial movies
bhvRange = [0 1];
frameRange = [1:152 457:760];
Widefield_SaveToAvi(newEyeTrace(:,:,frameRange), 'pupilExampleSmall', 15, 'gray', bhvRange, fPath, 4)

