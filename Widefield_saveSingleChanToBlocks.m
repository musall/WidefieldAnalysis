function Widefield_saveSingleChanToBlocks(cPath, recType, nrBlocks, overlap)
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
disp('Single-channel recording - No hemodynamic correction');
disp(datestr(now));

if strcmpi(recType, 'bpod')
    disp('BpodImager settings');
    opts.stimLine = 6; %analog line that contains stimulus trigger.
    opts.trigLine = 8; %analog lines for blue and violet light triggers.
elseif strcmpi(recType, 'widefield')
    disp('WidefieldImager settings');
    opts.stimLine = 3; %analog line that contains stimulus trigger.
    opts.trigLine = 6; %analog lines for blue and violet light triggers.
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
    blueData = Widefield_SplitChannels(opts,trials(1));
else
    blueData = Widefield_SplitVideo(opts,trials(1));
end
blueData = single(squeeze(blueData));
blueRef = fft2(median(blueData,3)); %blue reference for alignment
save([cPath 'blueRef.mat'],'blueRef');

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
blueFrameTimes = cell(1, fileCnt);
if ~exist([cPath 'blockData'], 'dir')
    mkdir([cPath 'blockData']);
end

if ~isempty(gpuDevice)
    blueRef = gpuArray(blueRef);
end

for iTrials = 1:fileCnt
    if opts.loadRaw
        [blueData,blueTimes] = Widefield_SplitChannels(opts,trials(iTrials));
    else
        [blueData,blueTimes] = Widefield_SplitVideo(opts,trials(iTrials));
    end
    
    if size(blueData,3) ~= size(hemoData,3)
        error(['Trial ' int2str(trials(iTrials)) ': Blue and hemo channels have uneven framecount'])
    end
    
    if ~isempty(gpuDevice)
        blueData = gpuArray(blueData);
    end
    
    %perform image alignment for both channels
    for iFrames = 1:size(blueData,3)
        [~, temp] = Widefield_dftregistration(blueRef, fft2(blueData(:, :, iFrames)), 10);
        blueData(:, :, iFrames) = abs(ifft2(temp));
    end
    blueData = gather(blueData);
    
    % keep avg for each trial to check if channels were separated correctly
    blueAvg(:,:,iTrials) = mean(blueData,3);
    
    %keep timestamps for all frames
    blueFrameTimes{iTrials} = blueTimes;

    if rem(iTrials,10) == 0
        fprintf(1, 'Loading session %d out of %d\n', iTrials,length(trials));
    end
    
    blueData = reshape(blueData, [], size(blueData,3));
    
    % save data in individual blocks. single file for each trial/block. Will delete those later.
    for iBlocks = 1:nrBlocks
        bBlock = blueData(blockInd{iBlocks}, :);
        save([cPath 'blockData' filesep 'blueBlock' num2str(iBlocks) '_Trial' num2str(iTrials)], 'bBlock', '-v6');
    end
end
clear blueData blueTimes blueRef
save([cPath 'trials.mat'],'trials'); %save trials so order of analysis is consistent

%save frametimes for blue/hemo trials
save([cPath 'blueFrameTimes.mat'],'blueFrameTimes', 'trials');

%save averages in case you need them later
save([cPath 'blueAvg.mat'],'blueAvg');

%take average over all trials for subsequent mean correction
blueAvg = mean(single(blueAvg),3);

%% subtract and divide each block by and compress with SVD
bU = cell(nrBlocks,1); bV = cell(nrBlocks,1);
for iBlocks = 1 : nrBlocks
    
    % rebuild current block from all trials
    Cnt = 0;
    allBlock = NaN(size(blockInd{iBlocks},1), size(cat(1,blueFrameTimes{:}),1), 'single');
    for iTrials = 1:fileCnt
        load([cPath 'blockData' filesep 'blueBlock' num2str(iBlocks) '_Trial' num2str(iTrials)], 'bBlock');
        allBlock(:, Cnt + (1:size(bBlock,2))) = single(bBlock);
        Cnt = Cnt + size(bBlock,2);
    end
    
    % compute dF/F   
    allBlock = bsxfun(@minus, allBlock, blueAvg(blockInd{iBlocks}));
    allBlock = bsxfun(@rdivide, allBlock, blueAvg(blockInd{iBlocks}));
    
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

% combine all dimensions and remove noise
[wbV, wbU] = Widefield_denoiseBlocks(bV, opts); 
save([cPath 'wbV.mat'],'wbU', 'wbV', 'blockInd', 'nrBlocks', '-v7.3'); %save whole-frame SVD results
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