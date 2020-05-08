function [blueData,blueTimes,hemoData,hemoTimes,stimOn,falseAlign] = Widefield_CheckChannels(opts,trialNr,blueHigh)

if exist('blueHigh','var')==0 || isempty(blueHigh)
    blueHigh = false; %usually blue frames have lower brightness
end
falseAlign = false;

%% load data and check channel identity
cFile = [opts.fPath filesep 'Analog_' num2str(trialNr) '.dat']; %current file to be read
[~,Analog] = Widefield_LoadData(cFile,'Analog'); %load analog data
Analog = double(Analog);

cFile = [opts.fPath filesep opts.fName '_' num2str(trialNr) '.dat']; %current file to be read
[header,data] = Widefield_LoadData(cFile,'Frames'); %load video data

%reshape data to compute mean frame intensities
dSize = size(data);
data = reshape(data,[],dSize(end));
temp = zscore(median(single(data)));
data = squeeze(reshape(data,dSize));
frameTimes = header(1:end-length(dSize)) * (86400*1e3); %extract frame times from header and convert to millisecond timestamps

if ~isnan(opts.trigLine)
    
    trace = Analog(opts.trigLine,:); %blue and violet light trigger channels
    trace = zscore(trace(1,end:-1:1) - trace(2,end:-1:1)); %invert and subtract to check color of last frame
    trace(round(diff(trace)) ~= 0) = 0; %don't use triggers that are only 1ms long
    lastBlue = find(trace > 1, 1);
    lastHemo = find(trace < -1,1);
    
    blueLast = lastBlue < lastHemo;
    if isempty(lastBlue) || isempty(lastHemo)
        warning(['Failed to find trigger signals. lastBlue: ' num2str(lastBlue) '; lastHemo: ' num2str(lastHemo) '; trialNr: ' num2str(trialNr)])
    end
    
    % find all triggers in stim line and choose the one that is on the longest as the true stimulus trigger 
    % apparently this line can be contaminated by noise
    stimOn = find(diff(double(Analog(opts.stimLine,:)) > 1500) == 1);
    ind = find((find(diff(double(Analog(opts.stimLine,:)) > 1500) == -1) - stimOn) > 2,1); %only use triggers that are more than 2ms long
    stimOn = stimOn(ind) + 1;

%     bFrame = find(temp < temp(end) + 6); %index for first black frame (assuming the last frame is really a dark frame)
    bFrame = find(temp < min(temp)*.75); %index for first black frame (assuming the last frame is really a dark frame)
    bFrame(bFrame < round(size(temp,2) / 2)) = []; %make sure black frame is in the second half of recording.
    bFrame = bFrame(1);
    blueInd = false(1,length(temp));

    if blueLast %last frame before black is blue
        if rem(bFrame,2) == 0 %blue frames (bFrame - 1) have uneven numbers
            blueInd(1:2:dSize(end)) = true;
        else
            blueInd(2:2:dSize(end)) = true;
        end
        lastFrame = size(Analog,2) - lastBlue; %index for end of last frame
    else %last frame before black is violet
        if rem(bFrame,2) == 0 %blue frames (bFrame - 2) have even numbers
            blueInd(2:2:dSize(end)) = true;
        else
            blueInd(1:2:dSize(end)) = true;
        end
        lastFrame = size(Analog,2) - lastHemo; %index for end of last frame
    end
    hemoInd = ~blueInd;
    hemoInd(1) = false; %make sure first frame is blue. That way we know that the first violet frame follows the first blue frame and the same index can be used for both.
    
    blueData = data(:,:,blueInd);
    hemoData = data(:,:,hemoInd);
    
    frameTimes = (frameTimes - frameTimes(bFrame - 1)) + lastFrame; %realign frameTime based on time of last non-dark frame
    blueTimes = frameTimes(blueInd);
    hemoTimes = frameTimes(hemoInd);

    if isempty(stimOn)
        blueStim = opts.preStim + 1 ; %if no stim was presented
        disp('Warning: No stimulus trigger found !')
        disp(['Current file: ' cFile])
    else
        blueStim = find((blueTimes - stimOn) > 0, 1);  %first blue frame after stimOn
    end
    selInd = blueStim - opts.preStim : blueStim + opts.postStim -1; %index for returned frames
            
    if min(selInd) < 1 || max(selInd) > length(blueTimes) || max(selInd) > length(hemoTimes)  %if reqested frames are not available in the current dataset
        selInd = 1 :  opts.preStim +  opts.postStim;
        if length(selInd) > length(blueTimes) || length(selInd) > length(hemoTimes) %if there are not enough frames available, just return maximum amount that is found in both channels
            selInd = 1:min([length(blueTimes) length(hemoTimes)]) ;
        end
        fprintf('Warning: Not enough pre- or postdata available. Returning first %d frames instead (%d requested).\n',length(selInd),opts.preStim +  opts.postStim)
        disp(['Current file: ' cFile])

        falseAlign = true; %return flag that alginment failed. Data should not be used in most cases.
    end
    
    blueTimes = blueTimes(selInd); %get blue frame times before and after stimOn
    blueData = blueData(:,:,selInd); %get blue frame data before and after stim on
    
    hemoTimes = hemoTimes(selInd); %get hemo frame times before and after stimOn
    hemoData = hemoData(:,:,selInd); %get hemo frame data before and after stim on
      
    %make sure each frame has a corresponding trigger signal
    sRate = floor(mean(diff(frameTimes)));
    frameJitter = sum(diff(frameTimes) > (sRate + 1));
    sRate = sRate - (1 - rem(sRate,2)); %make sure, inter-frame interval is uneven number
    blueInd = bsxfun(@minus, repmat(round(blueTimes),1,sRate),-floor(sRate/2):floor(sRate/2))'; %get index for each frameTime and get half the IFI before and after frame was acquired to check trigger signal.
    hemoInd = bsxfun(@minus, repmat(round(hemoTimes),1,sRate),-floor(sRate/2):floor(sRate/2))'; %get index for each frameTime and get half the IFI before and after frame was acquired to check trigger signal.
    
    blueTrig = zscore(Analog(opts.trigLine(1), blueInd(:))); %blue light trigger channel
    blueTrig = reshape(blueTrig, [], length(blueTimes));
    blueTrig = any(zscore(blueTrig) > 0.5);
    
    hemoTrig = zscore(Analog(opts.trigLine(2), hemoInd(:))); %blue light trigger channel
    hemoTrig = reshape(hemoTrig, [], length(hemoTimes));
    hemoTrig = any(zscore(hemoTrig) > 0.5);

    if sum(blueTrig) ~= length(blueTimes) || sum(hemoTrig) ~= length(hemoTimes)
        disp(['Potential frame time violations: ' num2str(frameJitter)])
        disp(['Confirmed blue trigs: ' num2str(sum(blueTrig)) '; blue frame times: ' num2str(length(blueTimes))])
        disp(['Confirmed hemo trigs: ' num2str(sum(hemoTrig)) '; hemo frame times: ' num2str(length(hemoTimes))])
        disp(['Current file: ' cFile])
    end

else %this finds the violet channel by its median brightness. blueHigh assumes that blue frames are brighter. Only use this in case of emergency !
    temp = diff(temp);
    temp = [-temp(1) temp];
    
    if blueHigh
        blueInd = temp > 0; %blue frames have higher brightness
    else
        blueInd = temp < 0; %blue frames have lower brightness
    end
    if round(mean(rem(find(blueInd),2))) == 1 %blue frames have uneven numbers
        blueInd = false(1,length(temp));
        blueInd(1:2:length(temp)) = true;
    else
        blueInd = false(1,length(temp)); %blue frames have uneven numbers
        blueInd(2:2:length(temp)) = true;
    end
    
    blueData = data(:,:,blueInd);
    blueData = blueData(:,:,1 : (opts.preStim + opts.postStim));
    blueTimes = frameTimes(blueInd);
    blueTimes = blueTimes(1 : (opts.preStim + opts.postStim));
    
    hemoData = data(:,:,~blueInd);
    hemoData = hemoData(:,:,1 : (opts.preStim + opts.postStim));
    hemoTimes = frameTimes(~blueInd);
    hemoTimes = hemoTimes(1 : (opts.preStim + opts.postStim));
    
end

%% plot result if requested
if opts.plotChans
    figure(50) %show result
    subplot(1,2,1); colormap gray
    imagesc(mean(blueData,3)); axis image
    title(['Blue frame average - Trial ' num2str(trialNr)])
    subplot(1,2,2);
    imagesc(mean(hemoData,3)); axis image
    title(['Hemo frame average - Trial ' num2str(trialNr)])
    drawnow;
end
end