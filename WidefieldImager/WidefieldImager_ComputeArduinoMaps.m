function WidefieldImager_ComputeArduinoMaps(Animal,Paradigm,Session,binSize,smth,rotateImage,nTrials,skipFirst)

%% Set basic variables
fclose('all');
% path = 'C:\data\WidefieldImager\Animals\';   %Widefield data path
% path = 'U:\space_managed_data\WidefieldImager\Animals\';   %Widefield data path
path = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\WidefieldImager\Animals\';   %Widefield data path

if ~exist('binSize','var')
    binSize = 4;                         %Size of box filter for spatial binning in pixels.
end

if ~exist('smth','var')
    smth = 2.5;                         %smoothing factor for gaussian smooth on maps
end

if ~exist('rotateImage','var')
    rotateImage = 0;                    %rotate image if the orientation is not as required
end

if ~exist('nTrials','var')
    nTrials = 0;                        %number of trials for phasemaps. If nCycles <=0 or > the number of cyles in the data. It will use all available cycles for a single phasemap.
end

if ~exist('settingInd','var')
    settingInd = [];
end

if ~exist('skipFirst','var')
    skipFirst = true;
end

opts.fName = 'Frames';
opts.fPath = [path Animal filesep Paradigm filesep Session];
opts.plotChans = true;
opts.trigLine = [6 7];
opts.blueHigh = true;
opts.stimLine = 3;
opts.preStim = 0.5;
opts.postStim = 5;
opts.alignRes = 10;
opts.alignRes = 10;

%% data source
Recs = dir([path Animal '\' Paradigm '\' Session]); %find recordings
disp(['Current path: ' path Animal '\' Paradigm '\' Session]); tic

%get overview of trials by looking at analog data
temp = ls([path Animal '\' Paradigm '\' Session '\Frames*mp4']);
temp1 = reshape(temp',1,numel(temp));
temp1 = strrep(temp1,'.mp4','    ');
temp = reshape(temp1,size(temp'))';clear temp1
Trials = sort(str2num(temp(:,8:end)));clear temp %identified trials
cFile = ls([path Animal '\' Paradigm '\' Session '\' Animal '*settings.mat']);

aStimData = [];
for iFiles = 1:size(cFile,1)
    load([path Animal '\' Paradigm '\' Session '\' cFile(iFiles,:)]);
    StimData.nTrials = str2double(StimData.handles.NrTrials);
    StimData.fCount = ones(1,StimData.nTrials).*iFiles;
    [StimData.VarNames,ind] = sort(StimData.VarNames);
    StimData.VarVals = StimData.VarVals(ind,:);
    aStimData = appendBehavior(aStimData,StimData); %stitch settings data together into one file
    clear StimData
end

if length(Trials) ~= sum(aStimData.nTrials) %compare recorded data to set amount of trials
    warning(['Recorded trials (' num2str(length(Trials)) ') are unequal to settings in visual stimulator (' num2str(aStimData.nTrials) ')'])
end

% get duration of individual trials. this is assumed to be constant by now.
stimType = aStimData.VarVals(ismember(aStimData.VarNames{1},'stimType'),1:length(Trials)); %Duration of a given trial
stimDur = aStimData.VarVals(ismember(aStimData.VarNames{1},'trialDuration'),1:length(Trials)); %Duration of a given trial
stimFreq = aStimData.VarVals(ismember(aStimData.VarNames{1},'cyclesPerSecond'),1:length(Trials)); %bar speed in a given trial
numCycles = stimDur(1:length(Trials))./(1./stimFreq(1:length(Trials))); %number of cycles in a given trial

allStimType = unique(stimType);
stimCnt = length(allStimType); %number of stimTypes in current recording.
for iCond = 1:stimCnt % check for number of cycles in each stimType and make sure they are uniform
    if length(unique(numCycles(stimType == allStimType(iCond)))) ~= 1
        error(['Number of cycles per trial is inconsistent for stimType ' num2str(allStimType(iCond)) '. This is not supported.'])
    end
end

nTrials(nTrials <= 0 | nTrials > length(Trials)/stimCnt) = length(Trials)/stimCnt; 
nTrials = unique(nTrials);
disp(['Trials per condition: ' num2str(length(Trials)/stimCnt) ' - Computing phasemaps every [' num2str(nTrials) '] trials']);
trialCnts = ones(stimCnt,length(nTrials));
avgData = cell(stimCnt,length(nTrials)); % averaged sequence for each condition and trialcount
avgTrace = cell(1,stimCnt); % averaged single cycle response for each condition and trialcount
fTransform = cell(stimCnt,length(nTrials)); % fourier transforms for each condition and trialcount 

%% get single trials and average to end up with one sequence
trialCnt = 0; %counter for trials (if 'Trials' is not continous)
condCnt = zeros(1,stimCnt); %counter for how many sequences were saved in each condition

for iTrials = Trials'
    %% load data
    
    [blueData,blueTimes,hemoData,hemoTimes,stimOn,falseAlign,sRate] = Widefield_SplitVideo(opts,iTrials,'mp4');
    
    Data = Widefield_HemoCorrect(blueData,hemoData,1:ceil(opts.preStim * sRate),5); %hemodynamic correction
    Data = squeeze(Data);

    cFile = [path Animal '\' Paradigm '\' Session '\Analog_' num2str(iTrials) '.dat']; %current file to be read
    [~,Analog] = Widefield_LoadData(cFile,'Analog'); %load analog data

    %% check for stimulus onset and compare measured stim frequency wih stimulator settings
    Analog = double(Analog);
    dSwitch = find(diff(Analog(3,:)) > 2000 == 1); %onset of arduino / PTB stimulation
    firstFrame = round(blueTimes(1)) + dSwitch; %timePoint of stimulus onset, based on trigger signal.
    
    mStimTrigs = length(dSwitch); %amount of produces stimuli
    mStimFreq = unique(round(1000 ./ diff(dSwitch),2)); %frequency of produced stimuli
    if length(mStimFreq) ~= 1
        error('Inconsistent timing between stimulus onset triggers. Check Analog signal and stimulator settings.')
    end
    
    %compare with stimulator settings
    if mStimTrigs ~= numCycles(trialCnt)
        warning(['Different number of cycles between triggers (' num2str(mStimTrigs) ') and stimulator settings (' num2str(numCycles(trialCnt)) '). Using trigger information for analysis.'])
    end
    if mStimFreq ~= stimFreq(iTrials)
        warning(['Different stimulus frequency between triggers (' num2str(mStimFreq) 'Hz) and stimulator settings (' num2str(stimFreq(iTrials)) 'Hz) . Using stimulator settings for analysis.'])
        mStimFreq = stimFreq(iTrials);
    end
    
    dSwitch = find(diff(Analog(6,:) > 1000) == 1); %onset of arduino / PTB stimulation

    sFrameDur = round(mean(diff(frameTimes - frameTimes(1)))); %average duration between acquired frames
    stimWin = 1000/mStimFreq/sFrameDur; %number of frames between two stim onsets
    dSwitch = diff(frameTimes - frameTimes(1)); %duration between all acquired frames
    if any(abs(dSwitch - sFrameDur) > sFrameDur*1.5)
        warning([num2str(sum(abs(dSwitch - sFrameDur) > sFrameDur*1.5)) ' lost camera frames in trial ' int2str(iTrials) '; something broken with camera settings?']);
        warning(['Average time difference between acquired frames: ' num2str(sFrameDur) ' ms']);
    end

    %% compute the amount of required frames and collect from data
    clear stimFrame
    for iStims = 1:length(firstFrame)
        temp = find(diff((frameTimes-firstFrame(iStims)) < 0));
        if ~isempty(temp)
            stimFrame(iStims) = temp; %index of frame that is closest to stimulus onset
        else
            stimFrame(iStims) = NaN;
        end
    end
    baseline = mean(Data(:,:,1:stimFrame(1)),3); %get average data before stim onset as baseline
    
    if skipFirst
        stimFrame(1) = [];  %remove first cycle because of onset transient
    end
    stimFrame = stimFrame - stimWin/2; %add half a cycle
    stimFrame(stimFrame < 1) = []; %reject too early stim onset
    mStimTrigs = length(stimFrame);
    
    Data = Data(:,:,round(stimFrame(1):1000/sFrameDur * (mStimTrigs/mStimFreq) + stimFrame(1) - 1)); %use data window at which stimulus was delievered
    stimFrame = stimFrame - stimFrame(1) + 1;
    
    binData = arrayResize(Data,binSize); %do spatial binning
    baseline = arrayResize(baseline,binSize); %do spatial binning
    clear Data
    
    %% running average
    for x = 1:length(nTrials)
        ind = stimType(iTrials) == allStimType;
        cTrialCnt = rem(condCnt(ind)+1,nTrials(x)); %current trial cycle for running average (reset to 1, when 'nTrials' is reached)
 
        if cTrialCnt == 1 || condCnt(ind) == 0 %start cycle for running average
            avgData{ind,x} = binData; %starting dataset for running average with set trialcount
        else
            avgData{ind,x} = (avgData{ind,x}.*cTrialCnt + binData) ./ (cTrialCnt+1); %produce running average
        end
        
        if cTrialCnt == 0 %reached requested trialcount. Compute fourier transform and increase counter
            temp = fft(avgData{ind,x},[],3);
            fTransform{ind,x}(trialCnts(ind,x),:,:) = squeeze(temp(:,:,mStimTrigs+1)); clear temp
            trialCnts(ind,x) = trialCnts(ind,x)+1; %increase counter for running average when required trialcount is reached
        end
    end
    
    binData = mean(reshape(binData,size(binData,1),size(binData,2),stimWin, length(stimFrame)),4); %average over single cycles
    
    if isempty(avgTrace{ind})
        avgTrace{ind} = binData;
        avgBaseline{ind} = baseline;
    else
        avgTrace{ind} = (avgTrace{ind}.*condCnt(ind) + binData) ./ (condCnt(ind)+1); %produce running average
        avgBaseline{ind} = (avgBaseline{ind}.*condCnt(ind) + baseline) ./ (condCnt(ind)+1); %produce running average for baseline
    end
    condCnt(ind) = condCnt(ind) +1; %increase condition counter
    disp(['Done loading trial ' int2str(iTrials) '/' int2str(max(Trials))]);
    clear binData
end

%% normalize avgTrace data
currentPath = [path Animal '\' Paradigm '\' Session];
save([currentPath filesep 'avgTrace.mat'],'avgTrace','avgBaseline');

for iCond = 1:stimCnt
    ind = stimType == allStimType(iCond);
    stimWin = mean(1000./stimFreq(ind)./sFrameDur); %number of frames between two stim onsets
    trace = mean(avgTrace{iCond}(:,:,stimWin/4:stimWin/2),3);
    for iFrames = 1:size(avgTrace{iCond},3)
        temp{iCond}(:,:,iFrames) = smooth2a((double(avgTrace{iCond}(:,:,iFrames)) - trace) ./ trace, smth);
    end
    
    trace = mean(temp{iCond}(:,:,stimWin/2:end),3);

    Widefield_SaveToAvi(temp{iCond}, ...
        ['20fps_StimType_' num2str(allStimType(iCond))],20,'jet',[min(min(trace)) max(max(trace))*2],currentPath,5)
    Widefield_SaveToAvi(temp{iCond}, ...
        ['5fps_StimType_' num2str(allStimType(iCond))],5,'jet',[min(min(trace)) max(max(trace))*2],currentPath,5)
end

%% get phase and magnitude maps and make figures
for iCond = 1:stimCnt
    for iTrials = 1:length(nTrials)
        for iRuns = 1:size(fTransform{iCond,iTrials},1)
            magMaps{iCond,iTrials}(:,:,iRuns) = imrotate(squeeze(abs(fTransform{iCond,iTrials}(iRuns,:,:))),rotateImage); %combined magnitude map.
            phaseMaps{iCond,iTrials}(:,:,iRuns) = imrotate(squeeze(angle(fTransform{iCond,iTrials}(iRuns,:,:))),rotateImage); %combined phase map
        end
        cMagMaps{iCond,iTrials} = squeeze(mean(magMaps{iCond,iTrials},3));
%       cMagMaps{iCond,iTrials} =(cMagMaps{iCond,iTrials}-min(cMagMaps{iCond,iTrials}(:)))./(max(cMagMaps{iCond,iTrials}(:))- min(cMagMaps{iCond,iTrials}(:))); %normalize between 0 and 1
        cPhaseMaps{iCond,iTrials} = squeeze(mean(phaseMaps{iCond,iTrials},3));
    end
end

%% save data and make figures
save([currentPath filesep 'fftData.mat'],'cMagMaps','cPhaseMaps','allStimType'); %save some data

set(0,'DefaultFigureWindowStyle','docked')
for iCond = 1:stimCnt
    for iTrials = 1:length(nTrials)
        figure('name',[Paradigm ' - ' Animal '; stimType = ' num2str(allStimType(iCond)) '; iTrials = ' num2str(iTrials)]);
        subplot(1,2,1)
        imagesc(smooth2a(cMagMaps{iCond,iTrials},smth,smth));axis square; colorbar; 
        caxis([0 700]); colormap jet; title([Paradigm ' - Normalized Magnitude']); 
%         freezeColors
        
        temp = cMagMaps{iCond,iTrials} > prctile(cMagMaps{iCond,iTrials}(:),0);
        subplot(1,2,2)
        vfsIm = imagesc(smooth2a(cPhaseMaps{iCond,iTrials},smth,smth)); colorbar; caxis([0 pi])
        set(vfsIm,'AlphaData',spatialFilterGaussian(temp,2)); axis square
        title([Paradigm ' - Thresholded Phase']);
%         colormap hsv; 
        export_fig([currentPath filesep 'StimType_' num2str(allStimType(iCond))],'-jpg')
        
    end
end
set(0,'DefaultFigureWindowStyle','normal')
