function WidefieldMapper_Analyze(Animal,Paradigm,varargin)
% Batch code to analyze data from the WidefieldMapper code. Saves down a 4D
% array 'binData' of x*y*Frames*trials for further analysis

%% Set basic variables
binSize = 4;                                %number of pixels to be integrated for first approximation of interesting areas
path = 'H:\WidefieldMapper\Animals\';  %Widefield data path
reAnalyze = true;                           %redo basic  analysis
stimOn = 21;                                %frame at which stimulus was presented
fclose('all');

%% data source
if ~isempty(varargin)
    Recs(1).name = varargin{1};
else
    Recs = dir([path Animal '\' Paradigm]); %find recordings
    Recs = Recs([Recs.isdir] & ~strncmpi('.', {Recs.name}, 1)); %remove folders that start with a dot
end

for iRecs = 1:length(Recs)
    temp = dir([path Animal '\' Paradigm '\' Recs(iRecs).name]);
    Cnt = 0; dCheck = false;

    for iFiles = 1:length(temp)
        if strcmpi(temp(iFiles).name,['binData_' int2str(binSize) '.mat'])
            Cnt = Cnt+1; %check for processed data
        end
        if strcmpi(temp(iFiles).name,'sTime_1.mat')
            dCheck = true; %check for imaging data
        end
    end
    clear temp
    
    if (Cnt ~= 1 && dCheck) || (reAnalyze  && dCheck)
        disp(['Current path: ' path Animal '\' Paradigm '\' Recs(iRecs).name]); tic
        temp = ls([path Animal '\' Paradigm '\' Recs(iRecs).name '\sTime*mat']); %get overview of trials by looking at analog data
        temp1 = reshape(temp',1,numel(temp));
        temp1 = strrep(temp1,'.mat','    ');
        temp = reshape(temp1,size(temp'))';clear temp1
        Trials = sort(str2num(temp(:,7:end)));clear temp %identified trials
        
        %% get some basic information on data structure
        load( [path Animal '\' Paradigm '\' Recs(iRecs).name '\Snapshot_1.mat']); %get snapshot to determine image size
        load( [path Animal '\' Paradigm '\' Recs(iRecs).name '\bTime_1.mat']); %get snapshot to determine image size
        load( [path Animal '\' Paradigm '\' Recs(iRecs).name '\sTime_1.mat']); %get snapshot to determine image size
        
        cFile = [path Animal '\' Paradigm '\' Recs(iRecs).name '\bFrames_1.dat']; %current file to be read
        [~,bData] = Widefield_LoadData_v2(cFile,'Frames'); %load video data
        [a,b] = size(Widefield_ShrinkImage_v2(squeeze(bData(:,:,1,1)),binSize)); %shrink one frame to make sure binData has the right size
        binData = zeros([a,b,length(bTime)+length(sTime),max(Trials)],'single'); clear a b bData
        
        %% get data, create a downsampled version and save in binData array
        for iTrials = Trials'
            cFile = [path Animal '\' Paradigm '\' Recs(iRecs).name '\bFrames_' num2str(iTrials) '.dat']; %current file to be read
            [~,bData] = Widefield_LoadData_v2(cFile,'Frames'); %load video data
            bData = squeeze(bData);
            
            cFile = [path Animal '\' Paradigm '\' Recs(iRecs).name '\sFrames_' num2str(iTrials) '.dat']; %current file to be read
            [~,sData] = Widefield_LoadData_v2(cFile,'Frames'); %load video data
            sData = squeeze(sData);
            
            load([path Animal '\' Paradigm '\' Recs(iRecs).name '\bTime_' num2str(iTrials) '.mat']); %load frametimes to compute frame rate
            temp = cat(1,bTime(:).AbsTime);temp = temp(:,end);

            %% basic analysis
            for iFrames = 1:size(bData,3)+size(sData,3)
                if iFrames > size(bData,3)
                    binData(:,:,iFrames,iTrials) = Widefield_ShrinkImage_v2(sData(:,:,iFrames-size(bData,3)),binSize); %get binned version
                else
                    binData(:,:,iFrames,iTrials) = Widefield_ShrinkImage_v2(bData(:,:,iFrames),binSize); %get binned version
                end
            end
            clear sData bData
        end
        
        %% save some data
        fPath = [path Animal '\' Paradigm '\' Recs(iRecs).name '\binData_' int2str(binSize) '.mat'];
        if ~isdir(['H' path(2:end) Animal '\' Paradigm '\' Recs(iRecs).name])
            mkdir(['H' path(2:end) Animal '\' Paradigm '\' Recs(iRecs).name]);
        end
        save(fPath,'binData','-v7.3'); %binned data, contains trials, X,Y and frames
        save(['H' fPath(2:end)],'binData'); %binned data, contains trials, X,Y and frames
        
        fPath = [path Animal '\' Paradigm '\' Recs(iRecs).name '\FullBaseFrames.mat'];
        baseFrames = squeeze(nanmean(nanmean(binData(:,:,1:stimOn,:),3),4));
        save(fPath,'baseFrames','-v7.3'); % baseline frame - can be used for normalization
        save(['H' fPath(2:end)],'baseFrames','-v7.3'); % baseline frame - can be used for normalization

        fPath = [path Animal '\' Paradigm '\' Recs(iRecs).name '\SecondHalfBaseFrames.mat'];
        baseFrames = squeeze(nanmean(nanmean(binData(:,:,floor(stimOn/2)+1:stimOn,:),3),4));
        save(fPath,'baseFrames','-v7.3'); % baseline frame with only second half of baseline - may help to reduce artefact from initial blue light response
        save(['H' fPath(2:end)],'baseFrames','-v7.3'); % baseline frame with only second half of baseline - may help to reduce artefact from initial blue light response

        fPath = [path Animal '\' Paradigm '\' Recs(iRecs).name '\AllFrames.mat'];
        baseFrames = squeeze(nanmean(nanmean(binData,3),4));
        save(fPath,'baseFrames','-v7.3'); % average over all frames
        save(['H' fPath(2:end)],'baseFrames','-v7.3'); % check if averaging over all frames may be helpful. Most likely it is not though.
        
        %% isolate different stim modes if present
        load([path Animal '\' Paradigm '\' Recs(iRecs).name '\StimModality.mat'])
        if length(unique(Stim)) > 1
            for iMods = unique(Stim)
                ind = Stim == iMods;
                
                fPath = [path Animal '\' Paradigm '\' Recs(iRecs).name '\pBinData_' int2str(binSize) '_' int2str(iMods) '.mat'];
                pBinData = binData(:,:,:,ind);
                save(fPath,'pBinData','-v7.3'); % baseline frame - can be used for normalization
                save(['H' fPath(2:end)],'pBinData','-v7.3'); % baseline frame - can be used for normalization
                
                disp(['Saving isolated modality data: pBinData_' int2str(binSize) '_' int2str(iMods) '.mat']);
                clear pBinData
            end
        end
        clear binData
        disp(['Finished path: ' path Animal '\' Paradigm '\' Recs(iRecs).name]);toc
        disp('=============================================')
    end
end