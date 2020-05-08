function Widefield_AlignFrames(path,fName,trigLine)

opts.fPath = 'H:\BpodImager\Animals\mSM23\SpatialDisc\07-Mar-2017';
opts.fName = 'Frames';
opts.trigLine = NaN;


%% load imaging data and do alignment
trials = Widefield_CheckFrameNrs(path,fName); %get trial numbers for all data files

cFile = [opts.fPath '\' opts.fName '_' num2str(trials(1)) '.dat']; %current file to be read
[~,data] = Widefield_LoadData(cFile,'Frames'); %load video data
refFrame = median(squeeze(data),3);
fftRef = fft2(refFrame);
disp(['Analyzing data from ' opts.fPath '; Loading "' opts.fName '_" files'])

for iTrials = 1:length(trials)
%% check for hemo channel if requested and plot result

if checkHemo

    [blueData,blueTimes,hemoData,hemoTimes] = Widefield_CheckChannels ...
    (fPath, fName, trialNr, trigLine, plotChans)


    %% split up into two data files and save again
    frameTimes = header(1:end-length(size(data))); %extract frame times from header
    data = squeeze(data);
    blueData = data(:,:,ind);
    blueTimes = frameTimes(ind);
    violetData = data(:,:,~ind);
    violetTimes = frameTimes(~ind);
    clear data
    
    figure(50) %show result
    subplot(1,2,1); colormap gray
    imagesc(mean(blueData,3)); axis image
    title(['Blue frame average - Trial ' num2str(trials(iTrials))])
    blueData = reshape(blueData,size(blueData,1),size(blueData,2),1,size(blueData,3));
    subplot(1,2,2);
    imagesc(mean(violetData,3)); axis image
    title(['Violet frame average - Trial ' num2str(trials(iTrials))])
    violetData = reshape(violetData,size(violetData,1),size(violetData,2),1,size(violetData,3));
    drawnow;
    
    
end
end




for iTrials = 1:length(trials)
    
    cFile = [path '\' fName '_' num2str(trials(iTrials)) '.dat']; %current file to be read
    [header,data] = Widefield_LoadData(cFile,'Frames'); %load video data
    data = gpuArray(data);
    
    temp = gpuArray(zeros(size(data),'uint16'));
    for iFrames = 1:size(data,4)
        [~, Greg] = Widefield_dftregistration(fftRef, fft2(data(:, :, 1, iFrames)));
        temp(:, :, 1, iFrames) = uint16(abs(ifft2(Greg)));
    end
    cData = gather(temp);
    



    fID = fopen([path '\alignFrames_' num2str(trials(iTrials)) '.dat'], 'wb'); %open binary stimulus file
    fwrite(fID,length(header),'double'); %write number of expected header values
    fwrite(fID,header(1:end-length(size(cData))),'double'); %write absolute timestamps of each frame
    fwrite(fID,size(cData),'double'); %write size of image data array
    fwrite(fID,cData,'uint16'); %write image data
    fclose(fID);
    
    if clearData %delete original data file to save local space
        delete([path '\' fName '_' num2str(trials(iTrials)) '.dat']);
    end
    disp(['Completed ' num2str(iTrials) '/' num2str(length(trials)) ' trials']) 
    
end
end