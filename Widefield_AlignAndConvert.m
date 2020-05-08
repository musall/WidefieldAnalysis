function Widefield_AlignAndConvert(fPath,fName,keepData)
% code to perform image alignment and hemodynamic correction for individual
% data files from WidefieldImager program.
%
% Usage: Widefield_AlignAndConvert(fPath,fName,keepData)
% Inputs: 
%           - fPath: path to data files
%           - fName: name to data files
%           - keepData: flag to indicate whether original data should be kept or deleted

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end

%% basic variables
opts.fPath = fPath; %construct path to data folder
opts.fName = fName; %name of files to be loaded. exptected structure is fName_trialNr
opts.plotChans = false; %plot results of blue/violet channel when loading raw data
opts.stimLine = 3; %line with stimulus trigger. This is used to align pre and poststim data
opts.trigLine = [6 7]; %trigger line for blue / violet LED respectively

opts.preStim = 1; %baseline
opts.postStim = 22; %poststim data
opts.alignRes = 10; %resolution for frame alignment. Default is 10, resulting in 1/10 of a pixel resolution.

%% do alignment / hemo correction
files = dir([opts.fPath opts.fName '*.dat']); %find files to convert

for iFiles = 1:length(files)
    
    [blueData,blueTimes,hemoData,~,~,~] = Widefield_SplitChannels(opts,iFiles);
    
    if size(blueData,3) ~= size(hemoData,3)
        error(['Trial ' int2str(trials(iTrials)) ': Blue and hemo channels have uneven framecount'])
    end
    
    if ispc
        blueData = gpuArray(single(blueData));
        hemoData = gpuArray(single(hemoData));
    else
        blueData = single(blueData);
        hemoData = single(hemoData);
    end
    
    blueRef = fft2(nanmean(blueData,3)); %blue reference for alignment
    hemoRef = fft2(nanmean(hemoData,3)); %blue reference for alignment
    
    for iFrames = 1:size(blueData,3) %perform image alignment for both channels
        [~, temp] = Widefield_dftregistration(blueRef, fft2(blueData(:, :, iFrames)), opts.alignRes);
        blueData(:, :, iFrames) = abs(ifft2(temp));
        
        [~, temp] = Widefield_dftregistration(hemoRef, fft2(hemoData(:, :, iFrames)), opts.alignRes);
        hemoData(:, :, iFrames) = abs(ifft2(temp));
    end
    
    [A,B,C] = size(blueData);
    blueDataAvg = mean(blueData(:,:,1:opts.preStim),3);
    blueData = bsxfun(@minus, blueData, blueDataAvg); % subtract baseline mean
    blueData = bsxfun(@rdivide, blueData, blueDataAvg); % divide by baseline mean
    blueData = reshape(blueData,[],C);
    
    hemoData = single(hemoData);
    hemoAvg = mean(hemoData(:,:,1:opts.preStim),3);
    hemoData = bsxfun(@minus, hemoData, hemoAvg); % subtract baseline mean
    hemoData = bsxfun(@rdivide, hemoData, hemoAvg); % divide by baseline mean
    hemoData = reshape(hemoData,[],C);
    
    %% smooth hemo data
    smoothFact = 5; % ensure kernel length is odd
    n = size(hemoData,2);
    cbegin = cumsum(hemoData(:,1:smoothFact-2),2);
    cbegin = bsxfun(@rdivide, cbegin(:,1:2:end), 1:2:(smoothFact-2));
    cend = cumsum(hemoData(:,n:-1:n-smoothFact+3),2);
    cend = bsxfun(@rdivide, cend(:,end:-2:1), (smoothFact-2:-2:1));
    
    hemoData = conv2(reshape(hemoData,[],C),ones(1,smoothFact)/smoothFact,'full'); %smooth trace with moving average of 'smoothFact' points
    hemoData = [cbegin hemoData(:,smoothFact:end-smoothFact+1) cend];
    
    %% perform regression and scale hemo to blue channel
    if ispc
        m = zeros(size(blueData,1),1,'gpuArray');
        b = zeros(size(blueData,1),1,'gpuArray');
    else
        m = zeros(size(blueData,1),1);
        b = zeros(size(blueData,1),1);
    end
    
    tic
    for iPix = 1:size(blueData,1)
        theta = [hemoData(iPix,:)' ones(size(blueData,2),1)] \ blueData(iPix,:)';
        m(iPix) = theta(1);
        b(iPix) = theta(2);
        
        if rem(iPix,10000) == 0
            fprintf(1, 'Current trial is %d out of %d\n', iPix,size(blueData,1));
            toc
        end
    end
    
    hemoData = bsxfun(@times, hemoData, m); % multiply hemoData with regression coefficient
    hemoData = bsxfun(@plus, hemoData, b); % add constant offset
    
    %% subtract hemo from blue, reshape and ensure there is no baseline offset
    blueData = blueData - hemoData; %subtract scaled hemoChannel from data
    blueData = reshape(blueData,A,B,C);
    blueDataAvg = median(blueData(:,:,1:opts.preStim),3);
    blueData = bsxfun(@minus, blueData, blueDataAvg); % correct baseline offset
    blueData = gather(blueData);
    blueData(blueData > 1) = 1; %clip data at too high values to ensure enough dynamic range with uint16s
    blueData(blueData < -1) = -1; %clip data at too low values to ensure enough dynamic range with uint16s
    
    blueData = mat2gray(blueData(:,:,2:end)); %normalize between 0 and 1
    blueData = uint16(blueData * 2^16); %convert to integers

    cFile = [opts.fPath 'ac' opts.fName '_' num2str(iFiles) '.dat'];
    Widefield_SaveData(cFile,blueData,blueTimes); %write new file for aligned blue data

    if ~keepData
        cFile = [opts.fPath opts.fName '_' num2str(iFiles) '.dat']; %current file to be deleted
        delete(cFile);
    end
end
    