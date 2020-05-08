opts.fPath = 'H:\BpodImager\Animals\mSM23\SpatialDisc\29-Nov-2016';
opts.fName = 'Frames';
opts.trigLine = NaN;
opts.checkHemo = true;
opts.plotChans = true;
opts.trigLine = NaN;
opts.SVDframeBin = 20;
opts.baselineFrames = 1:40;
opts.nSVD = 1000;

%% get single file information and allocate mov matrix
trials = Widefield_CheckFrameNrs(opts.fPath,opts.fName); %get trial numbers for all data files

for iTrials = 1:length(trials)
    cPath = [opts.fPath '\' opts.fName '_' num2str(trials(iTrials)) '.dat'];
    header = Widefield_LoadData(cPath,'Analog');
    frameCnt(iTrials) = header(end);
end

ind = frameCnt < (max(frameCnt) - max(frameCnt) /4); % only use trials that have at least 3/4 of all frames. Those should be the completed trials.
trials(ind) = []; frameCnt(ind) = [];

if max(header(end-3:end-2)) > 768 %resolution should not be equal or higher as 1024 in any dimension
    temp = floor(header(end-3:end-2) / ceil(max(header(end-3:end-2))/768)) ; % adjust movie resolution
else
    temp = header(end-3:end-2);
end
    
baseline = zeros(temp(1),temp(2), 'single');
mov = zeros(temp(1),temp(2),ceil(sum(frameCnt)/opts.SVDframeBin/2), 2,'single');

Cnt = zeros(1,2);
blueAvg = zeros(temp(1),temp(2), 'single');
hemoAvg = zeros(temp(1),temp(2), 'single');

%% perform image alignement for separate channels and collect data in mov matrix
for iTrials = 1:length(trials)
    
    cFile = [opts.fPath '\alignBlueFrames_' num2str(trials(iTrials)) '.dat'];
    if ~(exist(cFile,'file') == 2) %check if aligned data exists already
        [blueData,blueTimes,hemoData,hemoTimes] = Widefield_CheckChannels(opts,trials(iTrials));
        
        blueData = gpuArray(blueData);
        hemoData = gpuArray(hemoData);
        
        if iTrials == 1
            blueRef = fft2(median(blueData,3)); %blue reference for alignment
            hemoRef = fft2(median(hemoData,3)); %violet reference for alignment
        end
        
        for iFrames = 1:size(blueData,3)
            [~, temp] = Widefield_dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10);
            blueData(:, :, iFrames) = abs(ifft2(temp));
        end
        Widefield_SaveData(cFile,gather(blueData),blueTimes); %write new file for aligned blue channel.
        
        for iFrames = 1:size(hemoData,3)
            [~, temp] = Widefield_dftregistration(hemoRef, fft2(hemoData(:, :, iFrames)), 10);
            hemoData(:, :, iFrames) = abs(ifft2(temp));
        end
        Widefield_SaveData([opts.fPath '\alignHemoFrames_' num2str(trials(iTrials)) '.dat'],gather(hemoData),hemoTimes); %write new file for aligned hemo channel.
        
    else
        [header,blueData] = Widefield_LoadData(cFile,'Frames');
        blueTimes = header(1:end-length(size(blueData))) * (86400*1e3); %extract frame times from header and convert to millisecond timestamps
            
        [header,hemoData] = Widefield_LoadData([opts.fPath '\alignHemoFrames_' num2str(trials(iTrials)) '.dat'],'Frames');
        hemoTimes = header(1:end-length(size(hemoData))) * (86400*1e3); %extract frame times from header and convert to millisecond timestamps
    end
    
    if iTrials == 1
        mask = Widefield_getMask(squeeze(median(blueData,3)),25); %get dark threshold
    end
    
    blueAvg = (blueAvg*(iTrials-1) + mean(blueData(:,:,opts.baselineFrames),3))/iTrials; %create running average of baseline. Will be subtracted from mov matrix later.
    blueData(:,:,floor(size(blueData,3) / opts.SVDframeBin) * opts.SVDframeBin + 1:end) = []; %remove frames that are above divider
    blueData = reshape(blueData, size(blueData,1), size(blueData,2), [], opts.SVDframeBin); %reshape data and keep average in mov array
    mov(:,:,Cnt(1) + (1:size(blueData,3)),1) = gather(squeeze(mean(blueData,4)));
    Cnt(1) = Cnt(1) + size(blueData,3);
    
    hemoAvg = (hemoAvg*(iTrials-1) + mean(hemoData(:,:,opts.baselineFrames),3))/iTrials; %create running average of baseline. Will be subtracted from mov matrix later.
    hemoData(:,:,floor(size(hemoData,3) / opts.SVDframeBin) * opts.SVDframeBin + 1:end) = []; %remove frames that are above divider
    hemoData = reshape(hemoData, size(hemoData,1), size(hemoData,2), [], opts.SVDframeBin);
    mov(:,:,Cnt(2) + (1:size(hemoData,3)),2) = gather(squeeze(mean(hemoData,4)));
    Cnt(2) = Cnt(2) + size(hemoData,3);

    if rem(iTrials,10) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
end
blueAvg = gather(blueAvg);
hemoAvg = gather(hemoAvg);

save([opts.fPath '\blueAvg.mat'],'blueAvg');
save([opts.fPath '\hemoAvg.mat'],'hemoAvg');
clear blueData hemoData

%% subtract basline - do this in steps to reduce memory load
ind = 1:200:size(mov,3);
for x = 1:length(ind)
    if x == length(ind)
        mov(:,:,ind(x):end,1) = bsxfun(@minus, mov(:,:,ind(x):end,1), blueAvg); % subtract blue mean 
        mov(:,:,ind(x):end,1) = bsxfun(@rdivide, mov(:,:,ind(x):end,1), blueAvg); % divide by blue mean 
        mov(:,:,ind(x):end,2) = bsxfun(@minus, mov(:,:,ind(x):end,2), hemoAvg); % subtract violet mean 
        mov(:,:,ind(x):end,2) = bsxfun(@rdivide, mov(:,:,ind(x):end,2), hemoAvg); % divide by violet mean 
    else
        mov(:,:,ind(x):ind(x+1),1) = bsxfun(@minus, mov(:,:,ind(x):ind(x+1),1), blueAvg); % subtract blue mean 
        mov(:,:,ind(x):ind(x+1),1) = bsxfun(@rdivide, mov(:,:,ind(x):ind(x+1),1), blueAvg); % divide blue mean 
        mov(:,:,ind(x):ind(x+1),2) = bsxfun(@minus, mov(:,:,ind(x):ind(x+1),2), hemoAvg); % subtract violet mean 
        mov(:,:,ind(x):ind(x+1),2) = bsxfun(@rdivide, mov(:,:,ind(x):ind(x+1),2), hemoAvg); % subtract violet mean 
    end
end
mov = reshape(mov,size(mov,1)*size(mov,2),[]);
mov(mask(:),:) = 0;
mov = reshape(mov,size(baseline,1),size(baseline,2),[]);

%% compute svd
tic
disp('Computing SVD');
opts.nSVD = min(opts.nSVD, size(mov,3));
mov       = reshape(mov, [], size(mov,3));
COV       = mov' * mov/size(mov,1);
totalVar  = sum(diag(COV)); % total variance of data.
opts.nSVD = min(size(COV,1)-2, opts.nSVD);

[V, Sv]         = svd(gpuArray(double(COV)));
V = gather(V);
Sv = gather(Sv);

V               = V(:, 1:opts.nSVD);
Sv              = Sv(1:opts.nSVD, 1:opts.nSVD);
U               = single(normc(mov * V));
Sv              = single(diag(Sv));
clear COV mov
toc

size(U)

%% apply SVD to data
disp('Compute temporal component');
Cnt = zeros(1,2);

ind = frameCnt < (max(frameCnt) - max(frameCnt) /4); % only use trials that have at least 3/4 of all frames. Those should be complete trials.
trials(ind) = []; frameCnt(ind) = [];

blueV = zeros(opts.nSVD,sum(frameCnt)/2,'single');
hemoV = zeros(opts.nSVD,sum(frameCnt)/2,'single');
blueFirst = zeros(1,length(trials));

for iTrials = 1:length(trials)

    cFile = [opts.fPath '\alignBlueFrames_' num2str(trials(iTrials)) '.dat'];
    [header,blueData] = Widefield_LoadData(cFile,'Frames');
    blueTimes = header(1:end-length(size(blueData))) * (86400*1e3); %extract frame times from header and convert to millisecond timestamps
    [header,hemoData] = Widefield_LoadData([opts.fPath '\alignHemoFrames_' num2str(trials(iTrials)) '.dat'],'Frames');
    hemoTimes = header(1:end-length(size(hemoData))) * (86400*1e3); %extract frame times from header and convert to millisecond timestamps
  
    blueFirst(iTrials) = blueTimes(1) < hemoTimes(1); %check if the first frame in the recording is blue or violet
    
    if max(size(blueData(:,:,1))) > 768 %resolution should not be higher as 768 in any dimension
        blueData = arrayResize(blueData,ceil(max(size(blueData(:,:,1)))/768)) ; % adjust movie resolution
        hemoData = arrayResize(hemoData,ceil(max(size(hemoData(:,:,1)))/768)) ; % adjust movie resolution
    end
    blueData = single(blueData);    
    hemoData = single(hemoData);    
    
    blueData = bsxfun(@minus, blueData, blueAvg); % subtract mean
    blueData = bsxfun(@rdivide, blueData, blueAvg); % divide mean
    hemoData = bsxfun(@minus, hemoData, hemoAvg); % subtract mean
    hemoData = bsxfun(@rdivide, hemoData, hemoAvg); % divide mean

    blueV(:,Cnt(1) + (1:size(blueData,3))) = U' * reshape(blueData, [], size(blueData,3));
    hemoV(:,Cnt(2) + (1:size(hemoData,3))) = U' * reshape(hemoData, [], size(hemoData,3));
      
    Cnt(1) = Cnt(1) + size(blueData,3);
    Cnt(2) = Cnt(2) + size(hemoData,3);
    
    if rem(iTrials,10) == 0
        fprintf(1, 'Recompute session %d out of %d\n', iTrials,length(trials));
    end
    
end
U = reshape(U, size(blueAvg,1), size(blueAvg,2), []);

%% do hemodynamic correction
[correctV, T] = HemoCorrectLocal(U, SubSampleShift(blueV,1,2), hemoV, 40, [9 10], 5);

%% save output

%% re-create original trial structure and save
singleFrameCnt = unique(frameCnt)/2;

blueV = reshape(blueV,1000,singleFrameCnt,[]);
hemoV = reshape(hemoV,1000,singleFrameCnt,[]);
correctV = reshape(correctV,1000,singleFrameCnt,[]);

save([opts.fPath '\blueV.mat'],'blueV','U','Sv','totalVar');
save([opts.fPath '\hemoV.mat'],'hemoV','U','Sv','totalVar');
save([opts.fPath '\correctV.mat'],'correctV','U','T','Sv','totalVar');