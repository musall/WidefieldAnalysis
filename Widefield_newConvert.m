function Widefield_newConvert(cPath, Animal, Paradigm, Rec, saveRaw)

if ~exist('saveRaw','var') || isempty(saveRaw)
    saveRaw = false; %dont save raw data by default
end
    
if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

%construct path to data folder
opts.fPath = [cPath Animal filesep Paradigm filesep Rec];
try
    load([opts.fPath filesep 'handles.mat'],'handles'); %get handles from widefieldimager
    opts.frameRate = str2double(handles.FrameRate.String)/2; %frame rate of individual channels. With alternating illumination this is the absolute frameRate / 2.
catch
    opts.frameRate = 15;
end
opts.animal = Animal;
opts.paradigm = Paradigm;
opts.rec = Rec;
opts.fName = 'Frames'; %name of imaging data files.
opts.plotChans = false; %plot results for blue and violet channel when loading raw data. Will make the code break on the HPC.
opts.stimLine = 6; %analog line that contains stimulus trigger.
opts.barcodeLine = 7; %line for Bpod trialID.
opts.barcodeCheck = true; %check barcode signal in each trial to readout trialID from Bpod.
opts.trigLine = [8 9]; %analog lines for blue and violet light triggers.
opts.nSVD = 500; %number of dimensions from svd.
opts.maskThresh = 40; %threshold for data mask. This is in percent brightness. Higher values will create smaller masks.
opts.preStim = 3; %prestim duration in seconds.
opts.postStim = 4; %poststim duration in seconds.
opts.memLimit = 100; % memory limit for video data in workspace in gigabyte. Use frame subsampling to stay within limit.
opts.baselineFrames = 1:opts.frameRate; %1s baseline. this is used for dF/F analysis later.
opts.loadRaw = true; %flag to use raw data instead of .mj2 files.
opts.hemoSmooth = ceil(0.175 / (1/opts.frameRate)); %use 250ms window for moving average and smooth hemo data
opts.hemoSmooth = opts.hemoSmooth-1+mod(opts.hemoSmooth,2); % ensure kernel length is odd

%% check if .dat or .mj2 files are present. Use slower splitvideo code for latter case.
rawCheck = dir([opts.fPath filesep opts.fName '*dat']);
vidCheck = dir([opts.fPath filesep opts.fName '*mj2']);
analogCheck = dir([opts.fPath filesep 'Analog_*.dat']);

if size(rawCheck,1) == size(analogCheck,1)
elseif size(vidCheck,1) == size(analogCheck,1)
    opts.loadRaw = false;
    disp('Using .mj2 files instead of raw data.');
else
    error('Unequal number of imaging and analog data files. Aborted')
end

%% isolate good trials for analysis
tic
[trials, ~, trialCnt, bTrials] = Widefield_FindGoodTrials(opts);
fprintf('Using %d/%d recorded trials\n', length(trials), trialCnt);
disp(opts.fPath); disp(datestr(now));
frameCnt = repmat(round((opts.preStim + opts.postStim) * opts.frameRate), 1, length(trials));
maxFrameCnt = round((opts.preStim + opts.postStim) * opts.frameRate);
opts.falseAlign = false(1,length(trials));

cFile = [opts.fPath filesep 'mask.mat'];
if ~(exist(cFile,'file') == 2) %check if mask exists already
    if opts.loadRaw
        blueData = Widefield_SplitChannels(opts,trials(1));
    else
        blueData = Widefield_SplitVideo(opts,trials(1));
    end
    blueData = single(squeeze(blueData));
    trace = median(blueData,3); %smoothed mean image to create mask
    mask = Widefield_ManualMask(trace);
    trace(~mask) = 0;
    trace = smooth2a(double(trace),20,20); %smoothed mean image to create mask
    mask = ~imfill(trace > prctile(trace(:),opts.maskThresh),'holes');
    save([opts.fPath filesep 'mask.mat'],'mask');
else
    load([opts.fPath filesep 'mask.mat'],'mask');
end
clear trace header
    
if opts.loadRaw
    [blueData,~,hemoData] = Widefield_SplitChannels(opts,trials(1));
else
    [blueData,~,hemoData] = Widefield_SplitVideo(opts,trials(1));
end

% use indices to ensure data and mask have the same size
yRange = 1:size(mask,1); 
xRange = 1:size(mask,2);
blueData = blueData(yRange,xRange,:);
hemoData = hemoData(yRange,xRange,:);

blueData = single(squeeze(blueData));
blueRef = fft2(median(blueData,3)); %blue reference for alignment

hemoData = single(squeeze(hemoData));
hemoRef = fft2(median(hemoData,3)); %violet reference for alignment

blueAvg = zeros([size(mask), length(trials)],'uint16'); %average for mean correction. Collect single session averages to verify correct channel separation.
hemoAvg = zeros([size(mask), length(trials)],'uint16');

blueData = arrayShrink(blueData(:,:,1),mask,'merge'); %get one compressed frame for size estimate
info = whos('blueData');
exptSize = info.bytes * sum(frameCnt) * 2 / 2^30; %expected size of complete data set in gb.
frameFrac = ceil(exptSize / opts.memLimit); %only use fraction of frames to keep memory usage under control
fprintf(1,'MemLimit: %f gb ; Expected size: %f gb ; Selected fraction: %d; \n', opts.memLimit,exptSize,frameFrac);

mov = zeros(sum(~mask(:)),sum(floor(frameCnt/frameFrac)),2,'single'); %large mov matrix. Only use selected pixels later.
[~,selInd] = sort(rand(frameCnt(1),length(trials))); % get random index for each trial
selInd = selInd(1:floor(frameCnt(1)/frameFrac),:); % select first frameFrac numbers from index
save([opts.fPath filesep 'selInd.mat'],'selInd'); %save index for sub selection

info = whos('mov');
redSize = info.bytes / 2^30;
fprintf(1,'MemLimit: %f gb ; Expected size: %f gb ; Selected fraction: %d ; Reduced size: %f gb \n', opts.memLimit,exptSize,frameFrac,redSize);

%% perform image alignement for separate channels and collect data in mov matrix
Cnt = 0;
for iTrials = 1:length(trials)
    
    if opts.loadRaw
        [blueData,blueTimes,hemoData,hemoTimes,stimOn(iTrials),opts.falseAlign(iTrials)] = Widefield_SplitChannels(opts,trials(iTrials));
    else
        [blueData,blueTimes,hemoData,hemoTimes,stimOn(iTrials),opts.falseAlign(iTrials)] = Widefield_SplitVideo(opts,trials(iTrials));
    end
    blueData = blueData(yRange, xRange, :);
    hemoData = hemoData(yRange, xRange, :);
    
    if size(blueData,3) ~= size(hemoData,3)
        error(['Trial ' int2str(trials(iTrials)) ': Blue and hemo channels have uneven framecount'])
    end
    frameCnt(iTrials) = length(blueTimes);
    
    %perform image alignment for both channels
    for iFrames = 1:size(blueData,3)
        [~, temp] = Widefield_dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10);
        blueData(:, :, iFrames) = abs(ifft2(temp));
        
        [~, temp] = Widefield_dftregistration(hemoRef, fft2(hemoData(:, :, iFrames)), 10);
        hemoData(:, :, iFrames) = abs(ifft2(temp));
    end
    
    cFile = [opts.fPath filesep 'alignedBlueFrames_' num2str(trials(iTrials)) '.dat'];
    Widefield_SaveData(cFile,blueData,blueTimes); %write new file for aligned blue data
    cFile = [opts.fPath filesep 'alignedHemoFrames_' num2str(trials(iTrials)) '.dat'];
    Widefield_SaveData(cFile,hemoData,hemoTimes); %write new file for aligned hemo data
    
    % keep avg for each trial to check if channels were separated correctly
    blueAvg(:,:,iTrials) = median(blueData(:,:,opts.baselineFrames),3);
    hemoAvg(:,:,iTrials) = median(hemoData(:,:,opts.baselineFrames),3);
    
    %only use selected pixels from mask
    blueData = arrayShrink(blueData,mask,'merge'); 
    hemoData = arrayShrink(hemoData,mask,'merge');
    
    if frameFrac > 1 %subselect frames if required
        if opts.falseAlign(iTrials) %don't subselect if misaligned but cut down frames to match expected framecount
            if size(blueData,2) > size(selInd,1)
                blueData = blueData(:,1:size(selInd,1)); % use as many frames as you would have with subsampling
                hemoData = hemoData(:,1:size(selInd,1)); % use as many frames as you would have with subsampling
            end
        else
            blueData = blueData(:,selInd(:,iTrials)); % only use subselection of data
            hemoData = hemoData(:,selInd(:,iTrials)); % only use subselection of data
        end
    end

    mov(:,Cnt + (1:size(blueData,2)),1) = blueData;    
    mov(:,Cnt + (1:size(blueData,2)),2) = hemoData;
    Cnt = Cnt + size(blueData,2);

    if rem(iTrials,100) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
end
clear blueData hemoData blueRef hemoRef

%save all averages to check for channel separation
save([opts.fPath filesep 'allBlueAvg.mat'],'blueAvg');
save([opts.fPath filesep 'allHemoAvg.mat'],'hemoAvg');
 
%compute average across all sessions for later data correction
blueAvg = mean(single(blueAvg),3);
hemoAvg = mean(single(hemoAvg),3);
save([opts.fPath filesep 'blueAvg.mat'],'blueAvg');
save([opts.fPath filesep 'hemoAvg.mat'],'hemoAvg');

if size(mov,2) > Cnt(1)
    fprintf(1, 'Counter = %d; mov size = (%d,%d)\n', Cnt(1),size(mov,1),size(mov,2));
    fprintf(1, 'Removing %d/%d entries from mov matrix\n', size(mov,2) - Cnt(1),size(mov,2));
    mov(:,Cnt(1)+1:end,:) = []; %remove unused entries from mov matrix
end

%% subtract basline - do this in steps to reduce memory load
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

blueAvg = arrayShrink(blueAvg,mask,'split'); %recreate 2D image
hemoAvg = arrayShrink(hemoAvg,mask,'split'); %recreate 2D image
toc

%% compute svd
tic
disp('Computing SVD');
mov = reshape(mov,size(mov,1),[]); %merge channels
COV       = mov' * mov/size(mov,1); %compute covariance matrix over trials
totalVar  = sum(diag(COV)); % total variance of data.
opts.nSVD = min(size(COV,1)-2, opts.nSVD);

if size(COV,1) > 1e4
% if opts.nSVD < 1000 || size(COV,1) > 1e4
%     [V, Sv]  = eigs(double(COV), opts.nSVD, 'la');
    [V, Sv]  = rsvd(double(COV), opts.nSVD);
else
    [V, Sv]  = svd(COV);
    V        = single(V(:, 1:opts.nSVD));
    Sv       = single(Sv(1:opts.nSVD, 1:opts.nSVD));
end
clear COV

Sv = diag(Sv);
save([opts.fPath filesep 'firstV.mat'],'V','Sv','totalVar','-v7.3');
U = single(normc(mov * V));
clear V mov
toc; disp('done - loading data again to compute temporal component');

%% apply SVD to data
blueV = zeros(opts.nSVD,maxFrameCnt*length(trials),'single');
blueFrametimes = zeros(1,maxFrameCnt*length(trials));
hemoV = zeros(opts.nSVD,maxFrameCnt*length(trials),'single');
hemoFrameTimes = zeros(1,maxFrameCnt*length(trials));
Cnt = 0;

for iTrials = 1:length(trials)
    %load blue data and compute blueV
    cFile = [opts.fPath filesep 'alignedBlueFrames_' num2str(trials(iTrials)) '.dat'];
    [header,blueData] = Widefield_LoadData(cFile,'Frames');
    blueTimes = header(1:end-length(size(blueData))); %extract frame times from header and convert to milliseconds
    blueData = bsxfun(@minus, single(blueData), blueAvg); % subtract mean
    blueData = bsxfun(@rdivide, blueData, blueAvg); % divide mean
    blueData = arrayShrink(blueData,mask,'merge'); %only use selected pixels from mask
    delete(cFile); %delete aligned data file  
    
    %load hemo data and compute hemoV
    cFile = [opts.fPath filesep 'alignedHemoFrames_' num2str(trials(iTrials)) '.dat'];
    [header,hemoData] = Widefield_LoadData(cFile,'Frames');
    hemoTimes = header(1:end-length(size(hemoData))); %extract frame times from header and convert to milliseconds
    hemoData = bsxfun(@minus, single(hemoData), hemoAvg); % subtract mean
    hemoData = bsxfun(@rdivide, hemoData, hemoAvg); % divide mean
    hemoData = arrayShrink(hemoData,mask,'merge'); %only use selected pixels from mask
    delete(cFile); %delete aligned data file
    
    %check if trial has complete frame count and add NaNs to beginning
    preIdx = (opts.preStim * opts.frameRate + 1) - (find((blueTimes - stimOn(iTrials)) > 0, 1)); %check if sufficient baseline frames are present
    blueData = [NaN(size(blueData,1),preIdx) blueData];
    hemoData = [NaN(size(hemoData,1),preIdx) hemoData];
    blueTimes = [NaN(preIdx, 1); blueTimes];
    hemoTimes = [NaN(preIdx, 1); hemoTimes];
    
    %check frame count again and add NaNs to the end
    if size(blueData,2) < maxFrameCnt %if frames are missing at the end, add some NaNs
        postIdx = maxFrameCnt - size(blueData,2);
        blueData = [blueData NaN(size(blueData,1),postIdx)];
        hemoData = [hemoData NaN(size(hemoData,1),postIdx)];
        blueTimes = [blueTimes; NaN(postIdx, 1)];
        hemoTimes = [hemoTimes; NaN(postIdx, 1)];
    elseif size(blueData,2) > maxFrameCnt %if too many frames, limit to requested trial duration
        blueData = blueData(:, 1 : maxFrameCnt);
        hemoData = hemoData(:, 1 : maxFrameCnt);
        blueTimes = blueTimes(1 : maxFrameCnt);
        hemoTimes = hemoTimes(1 : maxFrameCnt);
    end
    
    % get blueV and hemoV
    blueV(:,Cnt + (1:size(blueData,2))) = U' * blueData;
    blueFrametimes(1,Cnt + (1:length(blueTimes))) = blueTimes;
    hemoV(:,Cnt + (1:size(hemoData,2))) = U' * hemoData;
    hemoFrameTimes(1,Cnt + (1:length(hemoTimes))) = hemoTimes;
    
    Cnt = Cnt + size(blueData,2);
    
    if rem(iTrials,100) == 0
        fprintf(1, 'Recompute session %d out of %d\n', iTrials,length(trials));
    end
end

%% re-create original matrix structures and save
U = arrayShrink(U,mask,'split'); %recreate spatial components
blueV = reshape(blueV,opts.nSVD,maxFrameCnt,[]);
blueFrametimes = reshape(blueFrametimes,maxFrameCnt,[]);
save([opts.fPath filesep 'blueV.mat'],'blueV','U','blueFrametimes','Sv','totalVar','trials','bTrials','-v7.3');

hemoV = reshape(hemoV,opts.nSVD,maxFrameCnt,[]);
hemoFrameTimes = reshape(hemoFrameTimes,maxFrameCnt,[]);
save([opts.fPath filesep 'hemoV.mat'],'hemoV','U','hemoFrameTimes','Sv','totalVar','trials','bTrials','-v7.3');

U = U(:, :, 1:200);
blueV = blueV(1:200, :, :);
hemoV = hemoV(1:200, :, :);
           
[Vc, regC, T, hemoVar] = Widefield_SvdHemoCorrect(U, blueV, hemoV, opts.frameRate);
save([opts.fPath filesep 'Vc.mat'],'Vc','U','blueFrametimes','Sv','totalVar','trials','bTrials','-v7.3');
save([opts.fPath filesep 'HemoCorrection.mat'],'regC','T', 'hemoVar')
save([opts.fPath filesep 'opts.mat'],'opts')
end