function Widefield_ConvertToSVD2(cPath, Animal, Paradigm, Rec)

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

opts.fPath = [cPath Animal filesep Paradigm filesep Rec]; %construct path to data folder
opts.animal = Animal;
opts.paradigm = Paradigm;
opts.rec = Rec;
opts.fName = 'Frames';
opts.trigLine = NaN;
opts.checkHemo = false;
opts.plotChans = false;
opts.stimLine = 6;
opts.barcodeLine = 7;
opts.barcodeCheck = true;
opts.trigLine = [8 9];
opts.nSVD = 2000;
opts.maskThresh = 35;
opts.preStim = 2;
opts.postStim = 2;
opts.memLimit = 150; % memory limit for video data in workspace in gigabyte. Use frame averaging to stay within limit.
opts.baselineFrames = 1:40;

%% isolate good trials for analysis
tic

[trials, ~, trialCnt] = Widefield_FindGoodTrials(opts);
fprintf('Using %d/%d recorded trials\n', length(trials), trialCnt);

[blueData,~,hemoData,~,~,~,opts.frameRate] = Widefield_SplitChannels(opts,trials(1));
opts.frameRate = round(1000/opts.frameRate); %frame rate of current dataset based on difference in timestamps in first selected trial
frameCnt = repmat(size(blueData,3) * 2, 1, length(trials));
falseAlign = false(1,length(trials));

cFile = [opts.fPath filesep 'mask.mat'];
if ~(exist(cFile,'file') == 2) %check if mask exists already
    blueData = Widefield_SplitChannels(opts,trials(1));
    blueData = single(squeeze(blueData));
    trace = median(blueData,3); %smoothed mean image to create mask
    mask = Widefield_ManualMask(trace);
    trace(~mask) = 0;
    trace = smooth2a(double(trace),20,20); %smoothed mean image to create mask
    mask = ~imfill(trace > prctile(trace(:),opts.maskThresh),'holes');
    blueHigh = false;
    save([opts.fPath filesep 'mask.mat'],'mask','blueHigh');
else
    load([opts.fPath filesep 'mask.mat']);
end
clear trace header
    
blueData = single(squeeze(blueData));
blueRef = fft2(median(blueData,3)); %blue reference for alignment

hemoData = single(squeeze(hemoData));
hemoRef = fft2(median(hemoData,3)); %violet reference for alignment

blueAvg = zeros([size(mask), length(trials)],'uint16'); %average for mean correction. Collect single session averages to verify correct channel separation.
hemoAvg = zeros([size(mask), length(trials)],'uint16');
    
blueData = arrayShrink(blueData(:,:,1),mask,'merge'); %get one compressed frame for size estimate
info = whos('blueData');
exptSize = info.bytes * sum(frameCnt) / 2^30; %expected size of complete data set in gb.
frameAvg = ceil(exptSize / opts.memLimit); %average across frames to keep memory usage under control
if frameAvg > 1 %adjust size estimate if averaging across frames
    frameAvg = find(rem(frameCnt(1)/2,frameAvg:frameCnt(1)) == 0,1) + frameAvg - 1; %make sure frameAvg can be divided by frameCnt/2
    exptSize = info.bytes * sum(frameCnt - (frameAvg*2)) / 2^30; % Subtract two 'frameAvg' from frameCnt because of random frame deletion later.
    mov = zeros(sum(~mask(:)),sum(frameCnt/2-frameAvg)/frameAvg,2,'single'); %large mov matrix. Only use selected pixels later.
    frameOffset = randi(frameAvg-1,1,length(trials));
    save([opts.fPath filesep 'frameOffset.mat'],'frameOffset');
else
    mov = zeros(sum(~mask(:)),sum(frameCnt/2),2,'single'); %large mov matrix. Only use selected pixels later.
end
fprintf(1,'MemLimit: %f ; Expected: %f ; FrameAvg: %d \n', opts.memLimit,exptSize,frameAvg);
Cnt = 0;

%% perform image alignement for separate channels and collect data in mov matrix
for iTrials = 1:length(trials)
    
    [blueData,blueTimes,hemoData,hemoTimes,~,falseAlign(iTrials)] = Widefield_SplitChannels(opts,trials(iTrials));
    if size(blueData,3) ~= size(hemoData,3)
        error(['Trial ' int2str(trials(iTrials)) ': Blue and hemo channels have uneven framecount'])
    end
    
    for iFrames = 1:size(blueData,3) %perform image alignment for both channels
        [~, temp] = Widefield_dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10);
        blueData(:, :, iFrames) = abs(ifft2(temp));
        
        [~, temp] = Widefield_dftregistration(hemoRef, fft2(hemoData(:, :, iFrames)), 10);
        hemoData(:, :, iFrames) = abs(ifft2(temp));
    end
    
    blueAvg(:,:,iTrials) = median(blueData(:,:,opts.baselineFrames),3);
    hemoAvg(:,:,iTrials) = median(hemoData(:,:,opts.baselineFrames),3);
    
    cFile = [opts.fPath filesep 'alignedBlueFrames_' num2str(trials(iTrials)) '.dat'];
    Widefield_SaveData(cFile,blueData,blueTimes); %write new file for aligned blue data
    cFile = [opts.fPath filesep 'alignedHemoFrames_' num2str(trials(iTrials)) '.dat'];
    Widefield_SaveData(cFile,hemoData,hemoTimes); %write new file for aligned hemo data
    
    blueData = arrayShrink(blueData,mask,'merge'); %only use selected pixels from mask
    hemoData = arrayShrink(hemoData,mask,'merge'); %only use selected pixels from mask
    
    if frameAvg > 1 %average frames if neccesarry.
        blueData(:,1:frameOffset(iTrials)) = []; %remove a random nr of frames to increase variability in the movie
        blueData(:,floor(size(blueData,2) / frameAvg) * frameAvg + 1:end) = []; %remove frames that are above average divider
        blueData = reshape(blueData, size(blueData,1), [], frameAvg); %reshape data and keep average in mov array
        blueData = squeeze(mean(blueData,3));
        
        hemoData(:,1:frameOffset(iTrials)) = []; %remove a random nr of frames to increase variability in the movie
        hemoData(:,floor(size(hemoData,2) / frameAvg) * frameAvg + 1:end) = []; %remove frames that are above average divider
        hemoData = reshape(hemoData, size(hemoData,1), [], frameAvg); %reshape data and keep average in mov array
        hemoData = squeeze(mean(hemoData,3));
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
mov = reshape(mov,size(mov,1),[]); %merge channels

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
    V        = single(V(:, 1:opts.nSVD));
    Sv       = single(Sv(1:opts.nSVD, 1:opts.nSVD));
end
clear COV

Sv        = diag(Sv);
toc;
disp('Saving V');
save([opts.fPath filesep 'firstV.mat'],'V','Sv','totalVar','-v7.3');
toc;

%% load mov again and compute U
if exptSize > 40 %had to delete 'mov' earlier to do svd computation. Reload.
    disp('done - loading mov matrix to compute U');
    if frameAvg > 1 %adjust size estimate if averaging across frames        
        mov = zeros(sum(~mask(:)),sum(frameCnt/2-frameAvg)/frameAvg,2,'single'); %large mov matrix. Using frame averaging.
    else
        mov = zeros(sum(~mask(:)),sum(frameCnt/2),2,'single'); %large mov matrix. No frame averaging.
    end
    Cnt = 0;
    
    for iTrials = 1:length(trials) %reload data, normalize and collect in mov array
        
        cFile = [opts.fPath filesep 'alignedBlueFrames_' num2str(trials(iTrials)) '.dat'];
        [~,blueData] = Widefield_LoadData(cFile,'Frames');
        blueData = bsxfun(@minus, single(blueData), blueAvg); % subtract mean
        blueData = bsxfun(@rdivide, blueData, blueAvg); % divide mean
        blueData = arrayShrink(blueData,mask,'merge'); %only use selected pixels from mask
        
        cFile = [opts.fPath filesep 'alignedHemoFrames_' num2str(trials(iTrials)) '.dat'];
        [~,hemoData] = Widefield_LoadData(cFile,'Frames');
        hemoData = bsxfun(@minus, single(hemoData), hemoAvg); % subtract mean
        hemoData = bsxfun(@rdivide, hemoData, hemoAvg); % divide mean
        hemoData = arrayShrink(hemoData,mask,'merge'); %only use selected pixels from mask
        
        if frameAvg > 1
            blueData(:,1:frameOffset(iTrials)) = []; %remove a random nr of frames to increase variability in the movie
            blueData(:,floor(size(blueData,2) / frameAvg) * frameAvg + 1:end) = []; %remove frames that are above average divider
            blueData = reshape(blueData, size(blueData,1), [], frameAvg); %reshape data and keep average in mov array
            
            hemoData(:,1:frameOffset(iTrials)) = []; %remove a random nr of frames to increase variability in the movie
            hemoData(:,floor(size(hemoData,2) / frameAvg) * frameAvg + 1:end) = []; %remove frames that are above average divider
            hemoData = reshape(hemoData, size(hemoData,1), [], frameAvg); %reshape data and keep average in mov array
        end
        
        mov(:,Cnt(1) + (1:size(blueData,2)),1) = squeeze(mean(blueData,3));
        mov(:,Cnt(1) + (1:size(blueData,2)),2) = squeeze(mean(hemoData,3));
        Cnt(1) = Cnt(1) + size(blueData,2);
        clear blueData hemoData
        
        if rem(iTrials,100) == 0
            fprintf(1, 'Loading aligned session %d out of %d\n', iTrials,length(trials));
        end
    end
end

if size(mov,2) > Cnt(1) 
    fprintf(1, 'Counter = %d; mov size = (%d,%d)\n', Cnt(1),size(mov,1),size(mov,2));
    fprintf(1, 'Removing %d/%d entries from mov matrix\n', size(mov,2) - Cnt(1),size(mov,2));
    mov(:,Cnt(1)+1:end,:) = []; %remove unused entries from mov matrix
end

mov = reshape(mov,size(mov,1),[]); %merge channels
U = single(normc(mov * V));
clear V mov
toc;
disp('done - loading data again to compute temporal component');

%% check if there was false alignment in any trial and reject those
if any(falseAlign)
    trials(falseAlign) = [];
    fprintf('Rejected %d/%d trials for too short baseline \n', sum(falseAlign), length(falseAlign))
end

%% apply SVD to data
blueV = zeros(opts.nSVD,median(frameCnt/2)*length(trials),'single');
blueFrametimes = zeros(1,median(frameCnt/2)*length(trials));
hemoV = zeros(opts.nSVD,median(frameCnt/2)*length(trials),'single');
hemoFrameTimes = zeros(1,median(frameCnt/2)*length(trials));
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

    blueV(:,Cnt + (1:size(blueData,2))) = U' * blueData;
    blueFrametimes(1,Cnt + (1:length(blueTimes))) = blueTimes;
    
    %load hemo data and compute hemoV
    cFile = [opts.fPath filesep 'alignedHemoFrames_' num2str(trials(iTrials)) '.dat'];
    [header,hemoData] = Widefield_LoadData(cFile,'Frames');
    hemoTimes = header(1:end-length(size(hemoData))); %extract frame times from header and convert to milliseconds
    hemoData = bsxfun(@minus, single(hemoData), hemoAvg); % subtract mean
    hemoData = bsxfun(@rdivide, hemoData, hemoAvg); % divide mean
    hemoData = arrayShrink(hemoData,mask,'merge'); %only use selected pixels from mask
    delete(cFile); %delete aligned data file
    
    hemoV(:,Cnt + (1:size(hemoData,2))) = U' * hemoData;
    hemoFrameTimes(1,Cnt + (1:length(hemoTimes))) = hemoTimes;
    
    Cnt = Cnt + size(blueData,2);
    
    if rem(iTrials,100) == 0
        fprintf(1, 'Recompute session %d out of %d\n', iTrials,length(trials));
    end
end

%% re-create original matrix structures and save
U = arrayShrink(U,mask,'split'); %recreate spatial components
blueV = reshape(blueV,opts.nSVD,median(frameCnt/2),[]);
blueFrametimes = reshape(blueFrametimes,median(frameCnt/2),[]);
save([opts.fPath filesep 'blueV.mat'],'blueV','U', 'blueFrametimes', 'Sv','totalVar','trials','-v7.3');

hemoV = reshape(hemoV,opts.nSVD,median(frameCnt/2),[]);
hemoFrameTimes = reshape(hemoFrameTimes,median(frameCnt/2),[]);
save([opts.fPath filesep 'hemoV.mat'],'hemoV','U', 'hemoFrameTimes', 'Sv','totalVar','trials','-v7.3');

[Vc, T] = Widefield_SvdHemoCorrect(U, blueV, hemoV, opts.frameRate);
save([opts.fPath filesep 'Vc.mat'],'Vc','U', 'blueFrametimes', 'Sv','totalVar','trials','-v7.3');
save([opts.fPath filesep 'hemoCorrection.mat'],'regC','T')

save([opts.fPath filesep 'opts.mat'],'opts') %save options that were used for conversion
end