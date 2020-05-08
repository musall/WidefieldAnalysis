function [U, Sv, V, hemoV, totalVar] = Widefield_ComputeSVD(fPath,fName,frameBin,baseLength)

opts.SVDframeBin = frameBin;
opts.dataPath = fPath;
opts.fileName = [fName '_'];
opts.baselineFrames = 1:baseLength;
% opts.nSVD = 500;
% opts.useGPU = false;
% opts.darkThresh = 25;
% opts.convertHemoChan = true;
% opts.hemoFileName = 'hemoFrames_';

%% get single file information and allocate mov matrix
trials = Widefield_CheckFrameNrs(opts.dataPath,opts.fileName); %get trial numbers for all data files
frameTimes = [];

for iTrials = 1:length(trials)
    cPath = [opts.dataPath '\' opts.fileName num2str(trials(iTrials)) '.dat'];
    header = Widefield_LoadData(cPath,'Analog');
    frameCnt(iTrials) = header(end);
    frameTimes = [frameTimes; [header(1:end-4), ones(1,size(header,1)-4)'*iTrials]];
end
save([fPath '\frameTimes.mat'],'frameTimes')

if max(header(end-3:end-2)) > 768 %resolution should not be equal or higher as 1024 in any dimension
    temp = floor(header(end-3:end-2) / ceil(max(header(end-3:end-2))/768)) ; % adjust movie resolution
else
    temp = header(end-3:end-2);
end
    
baseline = zeros(temp(1),temp(2), 'single');
mov = zeros(temp(1),temp(2),ceil(sum(frameCnt)/opts.SVDframeBin),'single');

%% load data and combine into large 'mov' matrix
disp('Load imaging data and subtract baseline');
tic
Cnt = 0;

% for iTrials = 1:length(trials)
for iTrials = 1
    
    cPath = [opts.dataPath '\' opts.fileName num2str(trials(iTrials)) '.dat'];
    [~,data] = Widefield_LoadData(cPath,'Frames'); data = squeeze(data);
    
    if max(header(end-3:end-2)) > 768 %resolution should not be higher as 768 in any dimension
        data = arrayResize(data,ceil(max(size(data(:,:,1)))/768)) ; % adjust movie resolution
    end
    data = single(data);
    data(:,:,floor(size(data,3) / opts.SVDframeBin) * opts.SVDframeBin + 1:end) = []; %remove frames that are above divider

    baseline = (baseline*(iTrials-1) + mean(data(:,:,opts.baselineFrames),3))/iTrials; %create running average of baseline. Will be subtracted from mov matrix later.
    
    data = reshape(data, size(data,1), size(data,2), [], opts.SVDframeBin);
    mov(:,:,Cnt + (1:size(data,3))) = squeeze(mean(data,4));
    Cnt = Cnt + size(data,3);
    
    if rem(iTrials,10) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
end


%% check for dark threshold
iThresh = opts.darkThresh;
Check = false;
h = figure('name','Check Threshold');
while ~Check
    imagesc(baseline);axis square; colormap gray; hold on
    contour(imfill(baseline > prctile(baseline(:),iThresh),'holes')); axis square; title(['Dark image index - Threshold: ' int2str(iThresh)])
    iThresh = prctile(baseline(:),iThresh); %threshold to detect darker part of the image. Pixels below threshold are excluded from further analysis.
    Wait = input('Happy with threshold? Enter "Y" or new threshold (0-100) to proceed \n','S');
    if strcmpi(Wait,'y')
        Check = true;
    elseif  ~isempty(str2num(Wait))
        Wait = str2num(Wait);
        disp(['Changed threshold to '  num2str(Wait)]);
        iThresh = Wait;
    end
end
close(h);
mask = imfill(baseline < iThresh, 'holes'); %mask to remove pixels

%% subtract basline - do this in steps to reduce memory load
ind = 1:200:size(mov,3);
for x = 1:length(ind)
    if x == length(ind)
        mov(:,:,ind(x):end) = bsxfun(@minus, mov(:,:,ind(x):end), baseline); % subtract mean here 
    else
        mov(:,:,ind(x):ind(x+1)) = bsxfun(@minus, mov(:,:,ind(x):ind(x+1)), baseline); % subtract mean here 
    end
end
mov = reshape(mov,size(mov,1)*size(mov,2),[]);
mov(mask(:),:) = 0;
mov = reshape(mov,size(baseline,1),size(baseline,2),[]);

%% compute svd
disp('Computing SVD');
tic
opts.nSVD = min(opts.nSVD, size(mov,3));
mov       = reshape(mov, [], size(mov,3));
COV       = mov' * mov/size(mov,1);
totalVar = sum(diag(COV)); % total variance of data.
opts.nSVD = min(size(COV,1)-2, opts.nSVD);

if opts.nSVD<1000 || size(COV,1)>1e4
    [V, Sv]          = eigs(double(COV), opts.nSVD);
else
    if opts.useGPU
        [V, Sv]         = svd(gpuArray(double(COV)));
        V = gather(V);
        Sv = gather(Sv);
    else
         [V, Sv]         = svd(COV);
    end
    V               = V(:, 1:opts.nSVD);
    Sv              = Sv(1:opts.nSVD, 1:opts.nSVD);
end

U               = single(normc(mov * V));
Sv              = single(diag(Sv));
clear COV mov
toc

%% apply SVD to data
disp('apply SVD to data');
Cnt = 0;
V = zeros(opts.nSVD,sum(frameCnt),'single');

for iTrials = 1:length(trials)
    
    cPath = [opts.dataPath '\' opts.fileName num2str(trials(iTrials)) '.dat'];
    [~,data] = Widefield_LoadData(cPath,'Frames'); data = squeeze(data);
    
    if max(size(data(:,:,1))) > 768 %resolution should not be higher as 768 in any dimension
        data = arrayResize(data,ceil(max(size(data(:,:,1)))/768)) ; % adjust movie resolution
    end
    data = single(data);
    
    data = bsxfun(@minus, data, baseline); % subtract mean
    V(:,Cnt + (1:size(data,3))) = U' * reshape(data, [], size(data,3));

    Cnt = Cnt + size(data,3);
    
    if rem(iTrials,10) == 0
        fprintf(1, 'Recompute session %d out of %d\n', iTrials,length(trials));
    end
end
U = reshape(U, temp(1),temp(2), []);
