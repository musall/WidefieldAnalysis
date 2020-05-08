tic
% opts.fPath = 'H:\BpodImager\Animals\mSM29\SpatialDisc\06-Apr-2017';
% opts.fPath = 'W:\data\BpodImager\Animals\mSM29\SpatialDisc\06-Apr-2017';
opts.fPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/mSM29/SpatialDisc/06-Apr-2017';
% opts.fPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/mSM27/SpatialDisc/06-Apr-2017';
opts.fName = 'Frames';
opts.trigLine = NaN;
opts.checkHemo = false;
opts.plotChans = false;
opts.trigLine = NaN;
opts.baselineFrames = 1:40;
opts.nSVD = 2000;
opts.maskThresh = 45;
opts.memLimit = 60; % memory limit for video data in workspace in gigabyte. Use frame averaging to stay within limit.

%% get single file information and allocate mov matrix
%get trial numbers for all data files
recs = dir([opts.fPath filesep opts.fName '*']);
for iRecs = 1:length(recs)
    a = textscan(recs(iRecs).name,'%s%f%s','Delimiter','_');
    trials(iRecs) = a{2};
end
trials = sort(trials);

for iTrials = 1:length(trials)
    cPath = [opts.fPath filesep opts.fName '_' num2str(trials(iTrials)) '.dat'];
    header = Widefield_LoadData(cPath,'Analog');
    frameCnt(iTrials) = header(end)/2; %only blue frames, so only use half the frameCnt
end

ind = frameCnt < (max(frameCnt)); % only use trials that have a complete framecount.
trials(ind) = []; frameCnt(ind) = [];

blueData = Widefield_CheckChannels(opts,1);
blueData = single(squeeze(blueData));
blueRef = fft2(median(blueData,3)); %blue reference for alignment
trace = smooth2a(median(double(blueData),3),10,10); %smoothed mean image to create mask
mask = ~imfill(trace > prctile(trace(:),opts.maskThresh),'holes');
save([opts.fPath filesep 'mask.mat'],'mask');
 
blueAvg = zeros(sum(~mask(:)),1, 'single'); %running average. Only use selected pixels later.
info = whos('blueAvg');
exptSize = info.bytes * sum(frameCnt) / 2^30; %expected size of complete data set in gb.
frameAvg = ceil(exptSize / opts.memLimit); %average across frames to keep memory usage under control
if frameAvg > 1 %adjust size estimate if averaging across frames
    exptSize = info.bytes * sum(frameCnt - frameAvg) / 2^30; % Subtract one 'frameAvg' from frameCnt because of random frame deletion later.
    mov = zeros(sum(~mask(:)),ceil(sum(frameCnt-frameAvg)/frameAvg),'single'); %large mov matrix. Only use selected pixels later.
    frameOffset = randi(frameAvg-1,1,length(trials));
else
    mov = zeros(sum(~mask(:)),ceil(sum(frameCnt)/frameAvg),'single'); %large mov matrix. Only use selected pixels later.
end
fprintf(1,'MemLimit: %f ; Expected: %f ; FrameAvg: %d \n', opts.memLimit,exptSize,frameAvg);
Cnt = zeros(1,2);

%% perform image alignement for separate channels and collect data in mov matrix
for iTrials = 1:length(trials)
    
    cFile = [opts.fPath filesep 'alignBlueFrames_' num2str(trials(iTrials)) '.dat'];
    if ~(exist(cFile,'file') == 2) %check if aligned data exists already
        [blueData,blueTimes] = Widefield_CheckChannels(opts,trials(iTrials));
        blueData = squeeze(blueData);
        
        for iFrames = 1:size(blueData,3)
            [~, temp] = Widefield_dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10);
            blueData(:, :, iFrames) = abs(ifft2(temp));
        end
        Widefield_SaveData(cFile,blueData,blueTimes); %write new file for aligned blue channel.
    else
        [~,blueData] = Widefield_LoadData(cFile,'Frames'); %load aligned data
    end
    
    blueData = arrayShrink(blueData,mask,'merge'); %only use selected pixels from mask
    blueAvg = (blueAvg*(iTrials-1) + mean(blueData(:,opts.baselineFrames),2))/iTrials; %create running average of baseline. Will be subtracted from mov matrix later.
    if frameAvg > 1
        blueData(:,1:frameOffset(iTrials)) = []; %remove a random nr of frames to increase variability in the movie
        blueData(:,floor(size(blueData,2) / frameAvg) * frameAvg + 1:end) = []; %remove frames that are above average divider
        blueData = reshape(blueData, size(blueData,1), [], frameAvg); %reshape data and keep average in mov array
    end
    
    mov(:,Cnt(1) + (1:size(blueData,2))) = squeeze(mean(blueData,3));
    Cnt(1) = Cnt(1) + size(blueData,2);

    if rem(iTrials,10) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
end
clear blueData

if size(mov,2) > Cnt(1) 
    fprintf(1, 'Counter = %d; mov size = (%d,%d)\n', Cnt(1),size(mov,1),size(mov,2));
    fprintf(1, 'Removing %d/%d entries from mov matrix\n', size(mov,2) - Cnt(1),size(mov,2));
    mov(:,Cnt(1)+1:end) = []; %remove unused entries from mov matrix
end

%% subtract basline - do this in steps to reduce memory load
ind = 1:200:size(mov,2);
for x = 1:length(ind)
    if x == length(ind)
        mov(:,ind(x):end,1) = bsxfun(@minus, mov(:,ind(x):end), blueAvg); % subtract blue mean 
        mov(:,ind(x):end,1) = bsxfun(@rdivide, mov(:,ind(x):end), blueAvg); % divide by blue mean 
    else
        mov(:,ind(x):ind(x+1),1) = bsxfun(@minus, mov(:,ind(x):ind(x+1)), blueAvg); % subtract blue mean 
        mov(:,ind(x):ind(x+1),1) = bsxfun(@rdivide, mov(:,ind(x):ind(x+1)), blueAvg); % divide blue mean 
    end
end
toc

blueAvg = arrayShrink(blueAvg,mask,'split'); %recreate 2D image
save([opts.fPath filesep 'blueAvg.mat'],'blueAvg');

%% compute svd
tic
disp('Computing SVD');
COV       = mov' * mov/size(mov,1);
totalVar  = sum(diag(COV)); % total variance of data.
opts.nSVD = min(size(COV,1)-2, opts.nSVD);

if exptSize > 40
    clear mov %clear mov matrix to give enough memory to compute svd
end

if opts.nSVD < 1000 || size(COV,1) > 1e4
    [V, Sv]  = eigs(double(COV), opts.nSVD, 'la');
else
    [V, Sv]  = svd(COV);
end

V         = V(:, 1:opts.nSVD);
Sv        = Sv(1:opts.nSVD, 1:opts.nSVD);
Sv        = single(diag(Sv));
clear COV
toc;

%% load mov again and compute U
recs = dir([opts.fPath filesep 'alignBlueFrames_*']);
for iRecs = 1:length(recs)
    a = textscan(recs(iRecs).name,'%s%f%s','Delimiter','_');
    trials(iRecs) = a{2};
end
trials = sort(trials);

if exptSize > 40 %had to delete 'mov' earlier to do svd computation. Reload.
    disp('done - loading mov matrix to compute U');
    Cnt = 0;
    
    for iTrials = 1:length(trials)
        
        cFile = [opts.fPath filesep 'alignBlueFrames_' num2str(trials(iTrials)) '.dat'];
        [header,blueData] = Widefield_LoadData(cFile,'Frames');
        blueData = single(blueData);
        
        blueData = bsxfun(@minus, blueData, blueAvg); % subtract mean
        blueData = bsxfun(@rdivide, blueData, blueAvg); % divide mean
        
        blueData = arrayShrink(blueData,mask,'merge'); %only use selected pixels from mask
        blueAvg = (blueAvg*(iTrials-1) + mean(blueData(:,opts.baselineFrames),2))/iTrials; %create running average of baseline. Will be subtracted from mov matrix later.
        if frameAvg > 1
            blueData(:,1:frameOffset(iTrials)) = []; %remove a random nr of frames to increase variability in the movie
            blueData(:,floor(size(blueData,2) / frameAvg) * frameAvg + 1:end) = []; %remove frames that are above average divider
            blueData = reshape(blueData, size(blueData,1), [], frameAvg); %reshape data and keep average in mov array
        end
        
        mov(:,Cnt(1) + (1:size(blueData,2))) = squeeze(mean(blueData,3));
        Cnt(1) = Cnt(1) + size(blueData,2);
        
        if rem(iTrials,10) == 0
            fprintf(1, 'Loading aligned session %d out of %d\n', iTrials,length(trials));
        end
    end
end

U = single(normc(mov * V));
clear V mov
toc;
disp('done - loading data again to compute temporal component');

%% apply SVD to data
blueV = zeros(opts.nSVD,max(frameCnt)*length(trials),'single');
blueFrametimes = zeros(1,max(frameCnt)*length(trials));
blueFirst = zeros(1,length(trials));
Cnt = 0;

for iTrials = 1:length(trials)

    cFile = [opts.fPath filesep 'alignBlueFrames_' num2str(trials(iTrials)) '.dat'];
    [header,blueData] = Widefield_LoadData(cFile,'Frames');
    blueTimes = header(1:end-length(size(blueData))) * (86400*1e3); %extract frame times from header and convert to milliseconds
    blueData = single(blueData);    
    
    blueData = bsxfun(@minus, blueData, blueAvg); % subtract mean
    blueData = bsxfun(@rdivide, blueData, blueAvg); % divide mean

    blueV(:,Cnt + (1:size(blueData,3))) = U' * reshape(blueData, [], size(blueData,3));
    blueFrametimes(1,Cnt + (1:length(blueTimes))) = blueTimes;
      
    Cnt = Cnt + size(blueData,3);
    
    if rem(iTrials,10) == 0
        fprintf(1, 'Recompute session %d out of %d\n', iTrials,length(trials));
    end
    
end

U = arrayShrink(U,mask,'split'); %recreate spatial components

%% re-create original trial structure and save
blueV = reshape(blueV,opts.nSVD,max(frameCnt),[]);
blueFrametimes = reshape(blueFrametimes,[],max(frameCnt));

save([opts.fPath filesep 'blueV.mat'],'blueV','U', 'blueFrametimes', 'Sv','totalVar','-v7.3');
