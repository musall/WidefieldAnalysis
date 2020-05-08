function Widefield_Preproc(Animal,Paradigm,Session)
% Batch code to analyze data from widefield imaging. This is created to
% analyse data from behaving animals in a bpod paradigm.

%% Set basic variables
fclose('all');
binSize = 4;                           %number of pixels for spatial binning to reduce data size
path = 'H:\BpodImager\Animals\';       %Widefield data path
reAnalyze = true;                      %flag to reAnalyze all data-sets. Otherwise only non-converted sets are analyzed.
stimOn = 40;                           %first frame after stimulus onset. everything before is baseline.
trimImage = true;                      %flag to indicate that dark parts of the image should be removed.
fileName = 'blueFrames_';

%% data source
if ~exist('Session','var') || isempty(Session)
    Recs = dir([path Animal '\' Paradigm]); %find recordings
else
    Recs(3).name = Session;
end

for iRecs = 3:length(Recs)
    temp = dir([path Animal '\' Paradigm '\' Recs(iRecs).name]);
    Cnt = 0; dCheck = false;

    for iFiles = 1:length(temp)
        if strcmpi(temp(iFiles).name,['binData_' int2str(binSize) '.mat'])
            Cnt = Cnt+1; %check for processed data
        end
        if strcmpi(temp(iFiles).name, [fileName '1.dat'])
            dCheck = true; %check for imaging data
        end
    end
    clear temp
    
    if (Cnt ~= 1 && dCheck) || (reAnalyze  && dCheck)
        disp(['Current path: ' path Animal '\' Paradigm '\' Recs(iRecs).name]); tic
        %% get overview of trials by looking at analog data
        temp = ls([path Animal '\' Paradigm '\' Recs(iRecs).name '\Analog*dat']); 
        temp1 = reshape(temp',1,numel(temp));
        temp1 = strrep(temp1,'.dat','    ');
        temp = reshape(temp1,size(temp'))';clear temp1
        Trials = sort(str2num(temp(:,8:end)))';clear temp %identified trials from file number
        
        %% load analog data and check barcodes to line up trial numbers
        for iTrials = 1:length(Trials)
            cFile = [path Animal '\' Paradigm '\' Recs(iRecs).name '\Analog_' num2str(Trials(iTrials)) '.dat']; %current file to be read
            [~,Analog] = Widefield_LoadData(cFile,'Analog'); %load analog data
            try
                bTrials(iTrials) = segmentVoltageAndReadBarcodes(double(Analog(7,:))./4095*5, 3, 6); %Bpod TrialNr, encoded from barcode signal
            catch
                bTrials(iTrials) = -1; % don't use this trial if barcode can't be recovered
            end
        end
        clear Analog
        if any(bTrials <= 0)
            disp([num2str(sum(bTrials<=0)) ' trials contain corrupted barcodes and are rejected from analysis'])
        end
        
        %% load behavior data and get indices for performed trials. This is to identify a completed trial to know about its framecount
        cFile = ls([path Animal '\' Paradigm '\' Recs(iRecs).name '\' Animal '_' Paradigm '*.mat']);
        load([path Animal '\' Paradigm '\' Recs(iRecs).name '\' cFile]);
        perfTrials = ismember(bTrials,find(~SessionData.DidNotLever));  %find performed trials - remove the rest from dataset        
         
        %% get some basic information on data structure
        load( [path Animal '\' Paradigm '\' Recs(iRecs).name '\Snapshot_1.mat']); %get snapshot to determine image size
        cFile = [path Animal '\' Paradigm '\' Recs(iRecs).name '\' fileName int2str(find(perfTrials,1)) '.dat']; %load first fully performed trial
        [~,Data] = Widefield_LoadData(cFile,'Frames'); %load video data
        disp(['Using ' num2str(sum(perfTrials)) '/'  num2str(length(perfTrials)) ' performed trials'])
            
        %% create mask to trim image
        if trimImage
            iThresh = 25; %default is 25% brightness as threshold
            temp  = squeeze(arrayResize(Data,binSize)); %get binned version
            temp = squeeze(mean(temp,3)); %mean image
            
            % check threshold for dark pixel rejection
            Check = false;
            h = figure('name','Check Threshold');
            while ~Check
                imagesc(temp);axis square; colormap gray; hold on
                contour(imfill(temp > prctile(temp(:),iThresh),'holes')); axis square; title(['Dark image index - Threshold: ' int2str(iThresh)])
                iThresh = prctile(temp(:),iThresh); %threshold to detect darker part of the image. Pixels below threshold are excluded from further analysis.
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
            mask = ~imfill(temp > iThresh, 'holes'); %mask to remove pixels
        else
            mask = false(size(arrayResize(squeeze(Data(:,:,1,1)),binSize))); %shrink one frame to make sure binData has the right size
        end
        clear temp
        binData = zeros([sum(sum(~mask)),size(Data,4),sum(perfTrials)],'uint16'); clear a b Data

        %% get data, create a downsampled version and save in binData array. Don't use incomplete trials.
        Cnt = 0;
        rCnt = 0;
        for iTrials = 1:length(bTrials)
            if perfTrials(iTrials) %check if trial was rejected, based on behavior or missing imaging data
                cFile = [path Animal '\' Paradigm '\' Recs(iRecs).name '\' fileName num2str(Trials(iTrials)) '.dat']; %current file to be read
                [~,Data] = Widefield_LoadData(cFile,'Frames'); %load video data
                Data = squeeze(Data);
                
                if size(Data,3) > 40
                    Cnt = Cnt+1;
                    Data  = uint16(arrayResize(Data,binSize)); %get binned version of frame stack
                    binData(:,:,Cnt) = arrayShrink(Data,mask); %merge dimensions and remove pixels based on mask
                    
                    %create running average
                    Data = double(Data);
                    baseline = mean(Data(:,:,1:stimOn),3);
                    for iFrames = 1:size(Data,3)
                        temp(:,:,iFrames) = (Data(:,:,iFrames) - baseline)./baseline; %normalized frames
                    end
                    if Cnt == 1
                        normTrialAvg = temp; %use these frames in the first trial
                    else
                        normTrialAvg = (normTrialAvg*(Cnt-1) + temp)/Cnt; %add frames to running average
                    end
                    clear temp
                else
                    rCnt = rCnt+1;
                    perfTrials(iTrials) = false;
                end
            end
            clear Data
        end
        disp(['Rejected ' num2str(rCnt) ' trials for containing 40 or less frames'])
        binData(:,:,Cnt+1:end) = []; %remove potential leftover entries in binData
        
        %% save some data
        fPath = [path Animal '\' Paradigm '\' Recs(iRecs).name '\binData_' int2str(binSize) '.mat'];
        if ~isdir([path Animal '\' Paradigm '\' Recs(iRecs).name])
            mkdir([path Animal '\' Paradigm '\' Recs(iRecs).name]);
        end
        save(fPath,'binData','mask','-v7.3'); %binned data, contains accepted pixels as X*Y x Frames x Trials. Mask is required to later reconstruct the actual frames using arrayShrink.
        
        fPath = [path Animal '\' Paradigm '\' Recs(iRecs).name '\normTrialAverage.mat'];
        save(fPath,'normTrialAvg','-v7.3'); % normalized average over all trials
        
        Widefield_SaveToAvi(normTrialAvg,'normTrialAvg_Video',5,'jet',[0 0.02],[path Animal '\' Paradigm '\'  Recs(iRecs).name])

        fPath = [path Animal '\' Paradigm '\' Recs(iRecs).name '\perfTrials.mat'];
        save(fPath,'perfTrials'); % index for performed trials
        
        clear binData baseFrames
        disp(['Finished path: ' path Animal '\' Paradigm '\' Recs(iRecs).name]);toc
        disp('=============================================')
    end
end
end