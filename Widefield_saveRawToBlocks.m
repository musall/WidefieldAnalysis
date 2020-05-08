function Widefield_saveRawToBlocks(cPath, recType, nrBlocks, overlap)
% code to collect all data from a given imaging experiment in individual
% blocks and perform dimensionality reduction.
% channels are either blue (1) or violet (2) which contains intrinsic
% signals. 
% cPath denotes the path to the data folder, containing raw video data and
% according timestamps. 
% nrBlocks is the total number of blocks used for the svd. Sqrt of blocks
% has to be an even number - nrBlocks is rounded down if need.
% overlap determines the number of pixels with which individual blocks are
% overlapping to avoid edge effects.

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

if ~exist('recType', 'var') || isempty(recType)
    recType = 'bpod';
end

if ~exist('nrBlocks', 'var') || isempty(nrBlocks)
    nrBlocks = 16;
end

% some variables for denoising
opts.blockDims = 200; %number of dimensions from SVD per block
opts.maxLag = 5; %lag for autocorrelation
opts.autoConfidence = 0.99; %confidence for autocorrelation test
opts.autoThresh = 1.5; %threshold for autocorrelation test
opts.snrThresh = 1.6; %threshold for SNR test
opts.svdMethod = 'randomized'; %method for blockwise SVD. either 'vanilla' or 'randomized'.

%construct path to data folder
opts.fPath = cPath; %build file path
opts.fName = 'Frames'; %name of imaging data files.
opts.plotChans = false; %plot results for blue and violet channel when loading raw data. Will make the code break on the HPC.
opts.loadRaw = true; %flag to use raw data instead of .mj2 files.
opts.verbosity = false; %flag to supress warnings from 'splitChannel / splitVideo' code.
opts.nrBlocks = nrBlocks; %keep number of blocks in structure for future reference
opts.overlap = overlap; %keep number overlap between blocks

disp('==============='); tic;
disp(opts.fPath); 
disp('Dual-channel recording - Apply hemodynamic correction');
disp(datestr(now));

if strcmpi(recType, 'bpod')
    disp('BpodImager settings');
    opts.stimLine = 6; %analog line that contains stimulus trigger.
    opts.trigLine = [8 9]; %analog lines for blue and violet light triggers.
elseif strcmpi(recType, 'widefield')
    disp('WidefieldImager settings');
    opts.stimLine = 3; %analog line that contains stimulus trigger.
    opts.trigLine = [6 7]; %analog lines for blue and violet light triggers.
end

if nrBlocks ~= floor(sqrt(nrBlocks))^2
    fprintf('Chosen nrBlocks (%d) cannot be squared. Using %d instead.\n', nrBlocks, floor(sqrt(nrBlocks))^2)
    nrBlocks = floor(sqrt(nrBlocks))^2;
end

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
save([cPath 'blueRef.mat'],'blueRef');
hemoData = single(squeeze(hemoData));
hemoRef = fft2(median(hemoData,3)); %violet reference for alignment
save([cPath 'hemoRef.mat'],'hemoRef');

%% get index for individual blocks
indImg = reshape(1:numel(blueRef),size(blueRef)); %this is an 'image' with the corresponding indices
blockSize = ceil((size(blueRef) + repmat(sqrt(nrBlocks) * overlap, 1, 2))/sqrt(nrBlocks)); %size of each block
blockInd = cell(1, nrBlocks);

Cnt = 0;
colSteps = (0 : blockSize(1) - overlap : size(blueRef,1)) + 1; %steps for columns
rowSteps = (0 : blockSize(2) - overlap : size(blueRef,2)) + 1; %steps for rows
for iRows = 1 : sqrt(nrBlocks)
    for iCols = 1 : sqrt(nrBlocks)
        
        Cnt = Cnt + 1;
        % get current block and save index as vector
        colInd = colSteps(iCols) : colSteps(iCols) + blockSize(1) - 1; 
        rowInd = rowSteps(iRows) : rowSteps(iRows) + blockSize(2) - 1;
        
        colInd(colInd > size(blueRef,1)) = [];
        rowInd(rowInd > size(blueRef,2)) = [];
        
        cBlock = indImg(colInd, rowInd);
        blockInd{Cnt} = cBlock(:);
        
    end
end
save([cPath 'blockInd.mat'],'blockInd');

%% perform image alignement for separate channels and collect data in mov matrix
blueAvg = zeros([size(blueData,1), size(blueData,2), fileCnt],'uint16'); %average for mean correction. Collect single session averages to verify correct channel separation.
hemoAvg = zeros([size(blueData,1), size(blueData,2), fileCnt],'uint16');
blueFrameTimes = cell(1, fileCnt);
hemoFrameTimes = cell(1, fileCnt);
if ~exist([cPath 'blockData'], 'dir')
    mkdir([cPath 'blockData']);
end

if ~isempty(gpuDevice)
    blueRef = gpuArray(blueRef);
    hemoRef = gpuArray(hemoRef);
end

for iTrials = 1:fileCnt
    if opts.loadRaw
        [blueData,blueTimes,hemoData,hemoTimes] = Widefield_SplitChannels(opts,trials(iTrials));
    else
        [blueData,blueTimes,hemoData,hemoTimes] = Widefield_SplitVideo(opts,trials(iTrials));
    end
    
    if size(blueData,3) ~= size(hemoData,3)
        error(['Trial ' int2str(trials(iTrials)) ': Blue and hemo channels have uneven framecount'])
    end
    
    if ~isempty(gpuDevice)
        blueData = gpuArray(blueData);
        hemoData = gpuArray(hemoData);
    end
    
    %perform image alignment for both channels
    for iFrames = 1:size(blueData,3)
        [~, temp] = Widefield_dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10);
        blueData(:, :, iFrames) = abs(ifft2(temp));
        
        [~, temp] = Widefield_dftregistration(hemoRef, fft2(hemoData(:, :, iFrames)), 10);
        hemoData(:, :, iFrames) = abs(ifft2(temp));
    end
    blueData = gather(blueData);
    hemoData = gather(hemoData);
    
    % keep avg for each trial to check if channels were separated correctly
    blueAvg(:,:,iTrials) = mean(blueData,3);
    hemoAvg(:,:,iTrials) = mean(blueData,3);
    
    %keep timestamps for all frames
    blueFrameTimes{iTrials} = blueTimes;
    hemoFrameTimes{iTrials} = hemoTimes;

    if rem(iTrials,10) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
    
    blueData = reshape(blueData, [], size(blueData,3));
    hemoData = reshape(hemoData, [], size(hemoData,3));
    
    % save data in individual blocks. single file for each trial/block. Will delete those later.
    for iBlocks = 1:nrBlocks
        bBlock = blueData(blockInd{iBlocks}, :);
        hBlock = hemoData(blockInd{iBlocks}, :);
        save([cPath 'blockData' filesep 'blueBlock' num2str(iBlocks) '_Trial' num2str(iTrials)], 'bBlock', '-v6');
        save([cPath 'blockData' filesep 'hemoBlock' num2str(iBlocks) '_Trial' num2str(iTrials)], 'hBlock', '-v6');
    end
end
clear blueData hemoData blueTimes hemoTiems blueRef hemoRef 
save([cPath 'trials.mat'],'trials'); %save trials so order of analysis is consistent

%save frametimes for blue/hemo trials
save([cPath 'blueFrameTimes.mat'],'blueFrameTimes', 'trials');
save([cPath 'hemoFrameTimes.mat'],'hemoFrameTimes', 'trials');

%save averages in case you need them later
save([cPath 'blueAvg.mat'],'blueAvg');
save([cPath 'hemoAvg.mat'],'hemoAvg');

%take average over all trials for subsequent mean correction
blueAvg = mean(single(blueAvg),3);
hemoAvg = mean(single(hemoAvg),3);
    
%% subtract and divide each block by and compress with SVD
bU = cell(nrBlocks,1); bV = cell(nrBlocks,1);
for iBlocks = 1 : nrBlocks
    
    % rebuild current block from all trials
    Cnt = 0;
    allBlock = NaN(size(blockInd{iBlocks},1), size(cat(1,blueFrameTimes{:}),1), 2, 'single');
    for iTrials = 1:fileCnt
        load([cPath 'blockData' filesep 'blueBlock' num2str(iBlocks) '_Trial' num2str(iTrials)], 'bBlock');
        load([cPath 'blockData' filesep 'hemoBlock' num2str(iBlocks) '_Trial' num2str(iTrials)], 'hBlock');
        
        allBlock(:, Cnt + (1:size(bBlock,2)), 1) = single(bBlock);
        allBlock(:, Cnt + (1:size(bBlock,2)), 2) = single(hBlock);
        Cnt = Cnt + size(bBlock,2);
    end
    
    % compute dF/F   
    allBlock(:,:,1) = bsxfun(@minus, allBlock(:,:,1), blueAvg(blockInd{iBlocks}));
    allBlock(:,:,1) = bsxfun(@rdivide, allBlock(:,:,1), blueAvg(blockInd{iBlocks}));
    allBlock(:,:,2) = bsxfun(@minus, allBlock(:,:,2), hemoAvg(blockInd{iBlocks}));
    allBlock(:,:,2) = bsxfun(@rdivide, allBlock(:,:,2), hemoAvg(blockInd{iBlocks}));
    
    % run randomized SVD on current block
    [bU{iBlocks}, bV{iBlocks}] = Widefield_compressSVD(reshape(allBlock,size(allBlock,1),[]), opts, false);
        
    if rem(iBlocks, round(nrBlocks / 5)) == 0
        fprintf(1, 'Loading block %d out of %d\n', iBlocks, nrBlocks);
    end
end

% save blockwise SVD data from both channels
save([cPath 'bV.mat'], 'bU', 'bV', 'blockInd', 'nrBlocks', '-v7.3');

% compute frame rate based on time difference between frames
frameDiff = diff(cat(1, blueFrameTimes{:}));
frameDiff(zscore(frameDiff) < -5) = [];
opts.sRate = round(1000 / mean(frameDiff));
save([cPath 'svdOpts.mat'], 'opts'); %save option struct

%% combine all dimensions and remove noise
[wbV, wbU] = Widefield_denoiseBlocks(bV, opts); 
save([cPath 'wbV.mat'],'wbU', 'wbV', 'blockInd', 'nrBlocks', '-v7.3'); %save whole-frame SVD results

% split channels and subtract means
wbV = reshape(wbV, size(wbV,1), [], 2);
wbV(:,:,1) = bsxfun(@minus, wbV(:,:,1), nanmean(wbV(:,:,1),2));
wbV(:,:,2) = bsxfun(@minus, wbV(:,:,2), nanmean(wbV(:,:,1),2));

% high-pass blueV and hemoV above 0.1Hz
[b, a] = butter(2,0.1/(opts.sRate), 'high');
wbV(:,:,1) = single(filtfilt(b,a,double(wbV(:,:,1)')))';
wbV(:,:,2) = single(filtfilt(b,a,double(wbV(:,:,2)')))';

% rebuild the denoised block SVD
wbV = reshape(wbV, size(wbV,1), []);
bV = (wbU * wbV);
bV = mat2cell(bV, ones(1, nrBlocks) * opts.blockDims, size(bV,2));

% run hemodynamic correction and save new Vc
[cV, regC, T, hemoExp] = Widefield_blockHemoCorrect(bV, bU, blueFrameTimes, opts);

% denoise corrected block data one last time.
[wcV, wcU] = Widefield_denoiseBlocks(cV, opts); 
save([cPath 'wcV.mat'], 'wcU', 'wcV', 'bU',  'blockInd', 'blueFrameTimes', 'nrBlocks', '-v7.3'); %save whole-frame corrected SVD results

fprintf('Blue rank : %d; Corrected rank: %d\n', size(wbV,1), size(wcV,1));
rmdir([cPath 'blockData' filesep], 's'); %remove temporary block data
toc;
disp('Done ===============');

end


%% reconstruct images based on individual blocks
% cV = (wcU * wcV);
% cV = mat2cell(cV, ones(1, nrBlocks) * opts.blockDims, size(cV,2));

% data1 = NaN(540*640,750);
% for iBlocks = 1 : size(cV,1)
%     temp = bU{iBlocks} * cV{iBlocks}(:, 1:750);
%     data1(blockInd{iBlocks},:) = nanmean(cat(3,data1(blockInd{iBlocks},:),temp),3);
% end
% data1 = reshape(data1,540,640,750);
% 
% newData = cat(4,data,data1);