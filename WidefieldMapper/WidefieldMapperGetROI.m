function mask = WidefieldMapperGetROI(Animal,Paradigm,Rec,snap,mask,rotateData)
% Code to analyze stacks that are created by WidefieldMapperAnalyze.

%% basic variables
binSize = 4;
path = 'H:\WidefieldMapper\Animals\'; %Widefield data path
smth = 5; %2-D smooth for plots
rThresh = 95; %treshold for ROI
Frames = 45; %frames to be used for analysis
cRange = []; %color range for selection image. Leave empty to scale.1
stimOn = 40; %frame at which stimulus was presented

%% load data
load([path Animal '\' Paradigm '\' Rec '\binData_' num2str(binSize) '.mat'])
% load([path Animal '\' Paradigm '\' Rec '\pBinData_' num2str(binSize) '_6.mat'])
% binData = pBinData; clear pBinData;

binData(isinf(binData)) = 0;
if ~exist('snap','var') || isempty(snap)
    snap = squeeze(mean(mean(binData(:,:,1:10,:),3),4));
end
if ~exist('rotateData','var') || isempty(rotateData)
    rotateData = 0;
elseif rem(rotateData,90) ~= 0
    rotateData = round(rotateData/90) * 90;
    disp(['RotateData is not a divider of 90. Rounding to ' num2str(rotateData) ' instead.']);
end
rotateData = rotateData/90; %make this a multiplier of 90 and use rot90 later.

%% correct for low frequency issues
for iTrials = 1:size(binData,4)
    Trace = squeeze(nanmean(nanmean(binData(:,:,:,iTrials),1),2));
    ind = Trace <= prctile(Trace,50); %define low-amp data as being below a certain percentile of all the data
    a = mean(Trace(ind)); b = std(Trace(ind));
    
    lInd = Trace >= a+b*2; %2STDs above low-amp data
    lDcross = find(diff(lInd) < 0); %instances where data was above low threshold and goes back to lower values
    lUcross = find(diff(lInd) > 0); %instances where data was above low threshold and goes back to lower values
    
    if lInd(1);lUcross = [1;lUcross];end
    if lInd(end);lDcross = [lDcross;length(lInd)];end
    
    
    hInd = Trace >= a+b*4; %instances where data passes 6STDs above low-amp data
    hCross = find(diff(hInd) < 0); %instances where data was above high threshold and goes back to lower values
    hCross = diff(hInd) == 1;
    rejFrames = [];
    
    for iThresh = 1:length(lUcross)
        
        if any(hInd(lUcross(iThresh):lDcross(iThresh)))
            rejFrames = [rejFrames lUcross(iThresh):lDcross(iThresh)];
        end
    end
    binData(:,:,rejFrames,iTrials) = NaN;
    
%     baseFrames = nanmean(binData(:,:,round(stimOn/2)+1:stimOn,iTrials),3);
    baseFrames = nanmean(binData(:,:,1:round(stimOn/2),iTrials),3);
    for iFrames = 1:size(binData,3)
        binData(:,:,iFrames,iTrials) = ((binData(:,:,iFrames,iTrials) - baseFrames) ./ baseFrames);
    end
end
binData = double(rot90(binData,rotateData));

%% Look at single frames and check if there is a sustained response
RejTrials = false;
if ~exist('mask','var') || isempty(mask)
    h1 = ShowStack(binData);
    
    if length(Frames) == 1
        Check = false;
        while ~Check
            Wait = input('Enter specific frame nr. ','S');
            if ~isnan(str2double(Wait))
                Frames = str2double(Wait);
                Check = true;
            end
        end
    end
    RejTrials = false;

    %% Draw ROI through average and check if user is happy
    [aInd,mask,h2] = FindROI_v4(binData,Frames(Frames>=stimOn),2,rThresh,smth,1,'k',cRange); %find ROI    
    Check = false;
    
    while ~Check
        Wait = input('Happy with ROI? (Y/N)  ','S');
        if strcmpi(Wait,'y')
            Check = true;
        elseif strcmpi(Wait,'n')
            close(h2);
            [aInd,mask,h2] = FindROI_v4(binData,Frames,2,rThresh,smth,1,'k',cRange); %find ROI
        elseif  ~isempty(str2num(Wait))
            close(h2);
            Wait = str2num(Wait);
            if length(Wait) == 1
                disp(['Changed threshold to '  num2str(Wait)]);
                rThresh = Wait; %treshold for ROI
            elseif length(Wait) == 2
                disp(['Changed caxis to ' num2str(Wait)]);
                cRange = Wait; %treshold for ROI
            end
            [aInd,mask,h2] = FindROI_v4(binData,Frames,2,rThresh,smth,1,'k',cRange); %find ROI
        end
    end
else
    h1 = figure;
    h2 = figure;
    aInd = contourc(double(mask),1);
    aInd(:,1) = [];
end
    
%% If result is insufficient, try to reject bad trials
RedoRejection = true;
tInd = true(1,size(binData,4)); %use all trials

if RejTrials
    if exist([path Animal '\' Paradigm '\' Rec '\RejTrials.mat'],'file')
        Check = false;
        while ~Check
            a = input('Reuse old rejection vector? (Y/N) ','S');
            if strcmpi(a,'y')
                load([path Animal '\' Paradigm '\' Rec '\RejTrials.mat'])
                Check = true; RedoRejection = false;
            elseif strcmpi(a,'n')
                Check = true; RedoRejection = true;
            end
        end
    end
    
    if RedoRejection
        close(h2);
        h2 = figure;colormap('jet')
        for x=1:size(binData,4)
            imagesc(squeeze(mean(binData(:,:,(Frames),x),3)));axis square;caxis([-0.015  0.015])
            Check = false;
            while ~Check
                UseTrial = input(['Use trial ' int2str(x) ' ? '],'S');
                if strcmpi(UseTrial,'Y')
                    UseTrial = true;
                    Check = true;
                elseif strcmpi(UseTrial,'N')
                    UseTrial = false;
                    Check = true;
                end
            end
            tInd(x) = UseTrial;
        end
    end
    close(h2);
    [aInd,mask,~] = FindROI_v4(binData(:,:,:,tInd),Frames,2,rThresh,smth,1,'k'); %find ROI
    save([path Animal '\' Paradigm '\' Rec '\RejTrials.mat'],'tInd')
end

%% get index for selected ROI and compute activity traces
xInd = [];
yInd = [];
xFind = find(sum(mask,1)>0);

for x = 1:length(xFind)
    yFind = find(mask(:,xFind(x)));
    for y = 1:length(yFind)
        xInd = [xInd xFind(x)];
        yInd = [yInd yFind(y)];
    end
end

roiTrace = nanmean(nanmean(Widefield_DimMerge_v1(binData,mask),1),3); %activity in region of interest
wholeTrace = squeeze(nanmean(nanmean(nanmean(binData,1),2),4)); %activity trace of all frames

%% load vessel image and add ROI to it
% figure;
% imagesc(snap); axis square; hold on; colormap gray
% plot(aInd(1,:)*binSize,aInd(2,:)*binSize,'c','linewidth',2)
% title([Animal ' - ' Rec])

h = figure('name',[Animal ' - ' Paradigm]);
subplot(2,2,1)
imagesc(rot90(snap,rotateData)); axis square; hold on; colormap(gca,'gray')
plot(aInd(1,:),aInd(2,:),'k','linewidth',2)
title([Animal ' - ' Rec])

subplot(2,2,2)
title([Animal ' - ' Rec])
plotMap = nanmean(nanmean(binData(:,:,Frames,tInd),3),4);
imagesc(smooth2a(plotMap,smth,smth));caxis([0 0.01]);
axis square;colormap(gca,'jet')
hold on;plot(aInd(1,:),aInd(2,:),'k','linewidth',2)
title([Animal ' - ' Rec])

subplot(2,2,3)
plot(wholeTrace);axis square;title('Whole image activity')
xlabel('# Frame');ylabel('delta R / R')

subplot(2,2,4)
plot(roiTrace);axis square;title('ROI activity')
xlabel('# Frame')

h.PaperUnits = 'inches';
set(h, 'PaperPosition', [0 0 15 15]);

% save figure and some data
savefig(h,[path Animal '\' Paradigm '\' Rec '\' Animal '_SelectedRegion.fig']);
saveas(h,[path Animal '\' Paradigm '\' Rec '\' Animal '_SelectedRegion.jpg'])
save([path Animal '\' Paradigm '\' Rec '\mask.mat'],'mask')
save([path Animal '\' Paradigm '\' Rec '\plotMap.mat'],'plotMap')

%% make movie and save
binData = (squeeze(nanmean(binData,4)));
for x = 1:size(binData,3)
    binData(:,:,x) = rot90(smooth2a(binData(:,:,x),smth,smth),2);
end
Widefield_SaveToAvi_v1(binData,[path Animal '\' Paradigm '\' Rec '_Video'],5,'jet',[0 0.015])

%%
close(h1)