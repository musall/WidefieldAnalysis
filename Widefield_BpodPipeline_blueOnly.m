opts.fPath = 'C:\data\WidefieldImager\Animals\mSM34\PhaseMap\24-Mar-2017';
% opts.fPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/mSM29/SpatialDisc/06-Apr-2017';
opts.fName = 'alignBlueFrames';
opts.trigLine = NaN;
opts.checkHemo = false;
opts.plotChans = true;
opts.trigLine = NaN;
opts.SVDframeBin = 1;
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

if max(header(end-2:end-1)) > 768 %resolution should not be equal or higher as 1024 in any dimension
    temp = floor(header(end-2:end-1) / ceil(max(header(end-3:end-2))/768)) ; % adjust movie resolution
else
    temp = header(end-2:end-1);
end
    
baseline = zeros(temp(1),temp(2), 'single');
mov = zeros(temp(1),temp(2),ceil((sum(frameCnt)/opts.SVDframeBin)), 1,'single');

Cnt = zeros(1,2);
blueAvg = zeros(temp(1),temp(2), 'single');

%% perform image alignement for separate channels and collect data in mov matrix
for iTrials = 1:length(trials)
    
    cFile = [opts.fPath '\alignBlueFrames_' num2str(trials(iTrials)) '.dat'];
    if ~(exist(cFile,'file') == 2) %check if aligned data exists already
        
        [header,blueData] = Widefield_LoadData([opts.fPath '\' opts.fName '_' num2str(trials(iTrials)) '.dat'],'Frames');
        blueTimes = header(1:end-length(size(blueData))); %extract frame times from header
        blueData = squeeze(blueData);
        
        if max(size(blueData(:,:,1))) > 768 %resolution should not be higher as 768 in any dimension
            blueData = arrayResize(blueData,ceil(max(size(blueData(:,:,1)))/768)) ; % adjust movie resolution
        end
        blueData = gpuArray(blueData);
        
        if iTrials == 1
            blueRef = fft2(median(blueData,3)); %blue reference for alignment
        end
        
        for iFrames = 1:size(blueData,3)
            [~, temp] = Widefield_dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10);
            blueData(:, :, iFrames) = abs(ifft2(temp));
        end
        Widefield_SaveData(cFile,gather(blueData),blueTimes); %write new file for aligned blue channel.
        
    else
        [header,blueData] = Widefield_LoadData(cFile,'Frames');
        blueTimes = header(1:end-length(size(blueData))); %extract frame times from header and convert to millisecond timestamps
    end
    
    if iTrials == 1
        mask = Widefield_getMask(squeeze(median(blueData,3)),25); %get dark threshold
    end
    
    blueAvg = (blueAvg*(iTrials-1) + mean(blueData(:,:,opts.baselineFrames),3))/iTrials; %create running average of baseline. Will be subtracted from mov matrix later.
    blueData(:,:,floor(size(blueData,3) / opts.SVDframeBin) * opts.SVDframeBin + 1:end) = []; %remove frames that are above divider
    blueData = reshape(blueData, size(blueData,1), size(blueData,2), [], opts.SVDframeBin); %reshape data and keep average in mov array
    mov(:,:,Cnt(1) + (1:size(blueData,3)),1) = gather(squeeze(mean(blueData,4)));
    Cnt(1) = Cnt(1) + size(blueData,3);

    if rem(iTrials,10) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
end
blueAvg = gather(blueAvg);

% save([opts.fPath '\blueAvg.mat'],'blueAvg');
clear blueData

%% subtract basline - do this in steps to reduce memory load
ind = 1:200:size(mov,3);
for x = 1:length(ind)
    if x == length(ind)
        mov(:,:,ind(x):end,1) = bsxfun(@minus, mov(:,:,ind(x):end,1), blueAvg); % subtract blue mean 
        mov(:,:,ind(x):end,1) = bsxfun(@rdivide, mov(:,:,ind(x):end,1), blueAvg); % divide by blue mean 
    else
        mov(:,:,ind(x):ind(x+1),1) = bsxfun(@minus, mov(:,:,ind(x):ind(x+1),1), blueAvg); % subtract blue mean 
        mov(:,:,ind(x):ind(x+1),1) = bsxfun(@rdivide, mov(:,:,ind(x):ind(x+1),1), blueAvg); % divide blue mean 
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

[V, Sv]         = svd((double(COV)));
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

blueV = zeros(opts.nSVD,sum(frameCnt),'single');
blueFrametimes = zeros(1,sum(frameCnt));
blueFirst = zeros(1,length(trials));

for iTrials = 1:length(trials)

    cFile = [opts.fPath '\alignBlueFrames_' num2str(trials(iTrials)) '.dat'];
    [header,blueData] = Widefield_LoadData(cFile,'Frames');
    blueTimes = header(1:end-length(size(blueData))) * (86400*1e3); %extract frame times from header and convert to milliseconds
    
    if max(size(blueData(:,:,1))) > 768 %resolution should not be higher as 768 in any dimension
        blueData = arrayResize(blueData,ceil(max(size(blueData(:,:,1)))/768)) ; % adjust movie resolution
    end
    blueData = single(blueData);    
    
    blueData = bsxfun(@minus, blueData, blueAvg); % subtract mean
    blueData = bsxfun(@rdivide, blueData, blueAvg); % divide mean

    blueV(:,Cnt(1) + (1:size(blueData,3))) = U' * reshape(blueData, [], size(blueData,3));
    blueFrametimes(1,Cnt(1) + (1:length(blueTimes))) = blueTimes;
      
    Cnt(1) = Cnt(1) + size(blueData,3);
    
    if rem(iTrials,10) == 0
        fprintf(1, 'Recompute session %d out of %d\n', iTrials,length(trials));
    end
    
end
U = reshape(U, size(blueAvg,1), size(blueAvg,2), []);

%% re-create original trial structure and save
singleFrameCnt = unique(frameCnt);
blueV = reshape(blueV,opts.nSVD,singleFrameCnt,[]);
blueFrametimes = reshape(blueFrametimes,[],singleFrameCnt);

save([opts.fPath '\blueV.mat'],'blueV','U', 'blueFrametimes', 'Sv','totalVar');
