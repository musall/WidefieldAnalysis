function Widefield_ComputeArduinoMaps(Animal,Session,winSize,smth,rotateImage,nTrials,pxPerMM)
%% for hindpaw example
% Widefield_ComputeArduinoMaps('mSM36','20-Oct-2017_1',4,2.5,0)
% path = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\WidefieldImager\Animals\';   %Widefield data path

%% Set basic variables
fclose('all');
% path = 'W:\data\WidefieldImager\Animals\';   %Widefield data path
% path = 'U:\space_managed_data\WidefieldImager\Animals\';   %Widefield data path
path = '\\churchlandNAS\homes\DOMAIN=CSHL\smusall\WidefieldImager\Animals\';   %Widefield data path
% path = 'H:\WidefieldImager\Animals\';   %Widefield data path
opts.fName = 'acFrames';
opts.fPath = [path Animal filesep 'MixedStim' filesep Session];
opts.plotChans = true;
opts.trigLine = [6 7];
opts.blueHigh = true;
opts.stimLine = 3;
opts.preStim = 0.5;
opts.alignRes = 10;
opts.alignRes = 10;
opts.hemoCorrect = false;
opts.sRate = 30;
smoothFact = 5;

if ~exist('winSize','var')
    winSize = 4;                        %size of the field in mm. This is to determine the spatial binning to get as close to 40pix/mm as possible. This is advised for the visual segmentation code.
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

if ~exist('pxPerMM','var')
    pxPerMM = [];                       %mapping from pixels to real space in pixels/mm. Typically, this is 165px/mm if fully zoomed in and 128px/mm if fully zoomed out. See the RulerPics in the WidefieldImager folder to check for yourself. 
end

%% data source
Recs = dir([path Animal '\MixedStim\' Session]); %find recordings
if isempty(Recs)
    path = 'C:\data\WidefieldImager\Animals\';   %alternate Widefield data path
    Recs = dir([path Animal '\MixedStim\' Session]); %find recordings
end
disp(['Current path: ' path Animal '\MixedStim\' Session]); tic

%get overview of trials by looking at analog data
useVideo = false;
temp = ls([path Animal '\MixedStim\' Session '\' opts.fName '*dat']);
temp1 = reshape(temp',1,numel(temp));
temp1 = strrep(temp1,'.dat','    ');
if isempty(temp1)
    temp = ls([path Animal '\MixedStim\' Session '\' opts.fName '*mj2']); 
    temp1 = reshape(temp',1,numel(temp));
    temp1 = strrep(temp1,'.mj2','    ');
    useVideo = true;
end
temp = reshape(temp1,size(temp'))';clear temp1
Trials = sort(str2num(temp(:,length(opts.fName)+2:end)));clear temp %identified trials

cFile = ls([path Animal '\MixedStim\' Session '\' Animal '*settings.mat']);
load([path Animal '\MixedStim\' Session '\' cFile]);

if length(Trials) ~= str2num(StimData.handles.NrTrials) %compare recorded data to set amount of trials
    trialOn = zeros(1,length(Trials));
    trialOff = zeros(1,length(Trials));
    
    for iTrials = 1:length(Trials)
        % Unequal trialcount. Search for missing trials in widefield data.
        cFile = [opts.fPath filesep opts.fName '_' num2str(iTrials) '.dat']; %current file to be read
        [header,~] = Widefield_LoadData(cFile,'Analog'); %load timestamps of imaging data
        
        trialOn(iTrials) = header(1) * 86400;
        trialOff(iTrials) = header(diff(header) < -1e3) * 86400;
    end
    trialOn(1) = [];
    trialOff(end) = [];

    trialDiff = trialOn - trialOff;
    idx = find((zscore(trialDiff) > 2)); % Find trials which have a too long ITI. Subsequent trial is most likely missed.
    StimData.TimeStamps(idx) = [];
    StimData.VarVals(:,idx) = []; 
end

if length(Trials) ~= (str2num(StimData.handles.NrTrials) - length(idx)) %compare recorded data to set amount of trials
    error(['Recorded trials (' num2str(length(Trials)) ') are unequal to settings in visual stimulator (' StimData.handles.NrTrials ')'])
else
    fprintf('Found %d trials that are missing from imaging data. Succesfully removed from StimData.\n',length(idx))
end

%get different bar direction and orentiations for individual trials
stimTypes = StimData.VarVals(strcmpi(StimData.VarNames,'stimType'),:); %get types of presented stimuli
allTypes = unique(stimTypes);

% get duration and bar speed for individual trials. this is assumed to be constant by now.
stimDur = StimData.VarVals(strcmpi(StimData.VarNames,'StimDuration') | strcmpi(StimData.VarNames,'trialDuration'),:); %Duration of a given trial
stimFreq = StimData.VarVals(strcmpi(StimData.VarNames,'cyclesPerSecond'),:); %stimulus frequency in a given trial
numCycles = unique(stimDur(Trials)./(1./stimFreq(Trials))); %number of cycles in current trial
opts.postStim = stimDur(1);

% check for number of cycles and trials in dataset
if length(numCycles) ~= 1
    error('Number of cycles per trial is inconsistent. This code is not meant to handle that.')
end
disp(['Trials per condition: ' num2str(length(Trials)/length(unique(stimTypes)))]);
avgData = cell(1,length(unique(stimTypes))); % averaged dataset for each condition
cData = cell(1,length(unique(stimTypes))); % current dataset for each condition

%% get single trials and average for each condition to end up with one sequence
tic;
condCnt = zeros(1,length(unique(stimTypes))); %counter for how many sequences were saved in each condition
blockCnt = 0; %counter for how many fft transforms were computed

for iTrials = Trials'
    %% load data    
    if strcmpi(opts.fName,'acFrames')
        cFile = [opts.fPath filesep opts.fName '_' num2str(iTrials) '.dat']; %current file to be read
        [header,Data] = Widefield_LoadData(cFile,'Frames'); %load video data
        frameTimes = header(1:find(diff(header) < -1e3)); %extract frame times from header and convert to millisecond timestamps
        sRate = 1000/round(median(diff(frameTimes))); %compute frame rate
    else
        if useVideo
            [blueData,~,hemoData,~,~,~,sRate] = Widefield_SplitVideo(opts,iTrials); %load mixed or blue only data from video files
        else
            [blueData,~,hemoData,~,~,~,sRate] = Widefield_SplitChannels(opts,iTrials); %load mixed or blue only data from .dat files
        end
        if opts.hemoCorrect
            Data = Widefield_HemoCorrect(blueData,hemoData,1:ceil(opts.preStim * sRate),5); %hemodynamic correction
        else
            Data = single(squeeze(blueData)); clear blueData
        end
    end
    Data = single(Data(:,:,1 : round((opts.preStim  + opts.postStim)* opts.sRate))); %trim data
    dataAvg = single(mean(Data(:,:,1:ceil(opts.preStim * sRate)),3));
    Data = (Data - dataAvg) ./ dataAvg;
         disp(iTrials);   
%     meanTrace(:,iTrials) = mean(reshape(Data,[],size(Data,3))); %get a mean average for the whole frame to check for suspicous activity

    cIdx = stimTypes(iTrials) == allTypes;
    condCnt(cIdx) = condCnt(cIdx) +1;
    cData{cIdx} = Data; %keep this in memory
        
    if condCnt(cIdx) == 1
        avgData{cIdx} = Data;
    else
        avgData{cIdx} = (avgData{cIdx}.*condCnt(cIdx) + Data) ./ condCnt(cIdx); %produce running average
    end
    
    %% perform fft after each block is over
    if ~any(condCnt < (blockCnt+1)) %if all conditions have a new recording, perform fft
        
        blockCnt = blockCnt + 1;

        temp = cat(4,cData{:});
        temp = reshape(temp(:,:,16:30,:),size(temp,1),size(temp,2),[]);
        temp = fft(temp,[],3); %combine data from all recordings and do fft
        fTransform{blockCnt,1} = squeeze(temp(:,:,2)); %find areas that respond to one stimulus
        fTransform{blockCnt,2} = squeeze(temp(:,:,length(cData) + 1)); %find areas that respond to all stimuli
        
        cIdx = ismember(allTypes, [5 6]);
        temp = cat(4,cData{cIdx});
        temp = reshape(temp(:,:,16:30,:),size(temp,1),size(temp,2),[]);
        temp = fft(temp,[],3); %combine data from audiovisual stimulation
        fTransform{blockCnt,3} = squeeze(temp(:,:,sum(cIdx) + 1)); %find areas that respond to both stimuli
        
        cIdx = allTypes > 6;
        temp = cat(4,cData{cIdx});
        temp = reshape(temp(:,:,16:30,:),size(temp,1),size(temp,2),[]);
        temp = fft(temp,[],3); %combine data from somatosensory stimulation
        fTransform{blockCnt,4} = squeeze(temp(:,:,sum(cIdx) + 1)); %find areas that respond to all stimuli
        
    end
end

%% temporal smooth for avg data
temp = cat(4,avgData{:});
temp = reshape(temp(:,:,1:30,:),size(temp,1),size(temp,2),[]);
[A,B,C] = size(temp);
temp = reshape(temp,[],C);

smoothFact = smoothFact-1+mod(smoothFact,2); % ensure kernel length is odd
n = size(temp,2);
cbegin = cumsum(temp(:,1:smoothFact-2),2);
cbegin = bsxfun(@rdivide, cbegin(:,1:2:end), 1:2:(smoothFact-2));
cend = cumsum(temp(:,n:-1:n-smoothFact+3),2);
cend = bsxfun(@rdivide, cend(:,end:-2:1), (smoothFact-2:-2:1));
temp = conv2(reshape(temp,[],C),ones(1,smoothFact)/smoothFact,'full'); %smooth trace with moving average of 'smoothFact' points
temp = [cbegin temp(:,smoothFact:end-smoothFact+1) cend];
temp = reshape(temp,A,B,[]);

%% %show average responses
labels = {'Visual' 'Forepaw(L)' 'Forepaw(R)' 'Hindpaw(L)' 'Hindpaw(R)' 'Trunk ' 'Snout(R)' 'Snout(L)'};
figure
for x = 1:8
    subplot(2,4,x)
    imagesc(arrayResize(mean(avgData{x+1}(:,:,18:20),3),1)); axis image
    title(labels{x}); colormap(inferno(256));
    caxis([0 0.3]);
end

%% show forepaw maps
figure
Cnt = 0;
for x = 5:-1:4
    Cnt = Cnt + 1;
    subplot(1,2,Cnt)
    imageScale(imrotate(smoothImg(mean(avgData{x+1}(:,:,15:25),3),1,5),-3));
    colormap(inferno(256)); caxis([0 0.3]); title([Animal ' - ' Rec ' - ' labels{x}]);
end