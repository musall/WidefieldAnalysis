function WidefieldSVD_ComputePhaseMaps(Animal,Session,winSize,smth,rotateImage,nTrials,pxPerMM)

%% Set basic variables
fclose('all');
path = 'H:\WidefieldImager\Animals\';   %Widefield data path
phaseMapSmth = 1;                       % smoothing factor for phaseMap. This is the first smoothing step to get more robust gradients. Too high values cause false 'areas'.
fName = 'blueV_5avg';

if ~exist('winSize','var')
    winSize = 4;                         %size of the field in mm. This is to determine the spatial binning to get as close to 40pix/mm as possible. This is advised for the visual segmentation code.
end

if ~exist('smth','var')
    smth = 1.5;                         %smoothing factor for gaussian smooth on maps
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
Recs = dir([path Animal '\PhaseMap\' Session]); %find recordings
if isempty(Recs); %sometimes recordings are stored on different hdd
    path = 'C:\data\WidefieldImager\Animals\';   %alternate Widefield data path
    Recs = dir([path Animal '\PhaseMap\' Session]); %find recordings
end
disp(['Current path: ' path Animal '\PhaseMap\' Session]); tic

%Load svd-compressed data. Contains all imaging data and timestamps in ms for each frame.
load([path Animal '\PhaseMap\' Session '\' fName '.mat']);

%load settings
cFile = ls([path Animal '\PhaseMap\' Session '\' Animal '*settings.mat']);
load([path Animal '\PhaseMap\' Session '\' cFile]);
Trials = 1 : str2num(StimData.handles.NrTrials);

%% get different bar direction and orentiations for individual trials
BarOrient = StimData.VarVals(strcmpi(StimData.VarNames,'BarOrient'),:); %get bar orientation for all trials
BarDirection = StimData.VarVals(strcmpi(StimData.VarNames,'BarDirection'),:); %get bar direction for all trials

iBarConds(1,:) = BarOrient == 1  & BarDirection == 0; % horizontal moving down
iBarConds(2,:) = BarOrient == 1  & BarDirection == 1; % horizontal moving up
iBarConds(3,:) = BarOrient == 0  & BarDirection == 0; % vertial, moving left to right
iBarConds(4,:) = BarOrient == 0  & BarDirection == 1; % vertical, moving right to left

% get duration and bar speed for individual trials. this is assumed to be constant by now.
StimDur = StimData.VarVals(strcmpi(StimData.VarNames,'StimDuration') | strcmpi(StimData.VarNames,'trialDuration'),:); %Duration of a given trial
barFreq = StimData.VarVals(strcmpi(StimData.VarNames,'cyclesPerSecond'),:); %bar speed in a given trial
numCycles = unique(StimDur(Trials)./(1./barFreq(Trials))); %number of cycles in current trial

% check for number of cycles and trials in dataset
if length(numCycles) ~= 1
    error('Number of cycles per trial is inconsistent. This code is not meant to handle that.')
end
nTrials(nTrials <= 0 | nTrials > length(Trials)/4) = length(Trials)/4; 
nTrials = [nTrials length(Trials)/4]; nTrials = unique(nTrials);
disp(['Trials per condition: ' num2str(length(Trials)/4) ' - Computing phasemaps from [' num2str(nTrials) '] trials']);
TrialCnts = ones(4,length(nTrials));
avgData = cell(4,length(nTrials)); % averaged sequence for each condition and trialcount
fTransform = cell(4,length(nTrials)); % fourier transforms for each condition and trialcount 

%% get single trials and average for each condition to end up with one sequence
condCnt = zeros(1,4); %counter for how many sequences were saved in each condition
for iTrials = Trials
    %% load data
    Data = svdFrameReconstruct(U, blueV(:,:,iTrials));
    frameTimes = blueFrametimes(:,iTrials);
    
    cFile = [path Animal '\PhaseMap\' Session '\Analog_' num2str(iTrials) '.dat']; %current file to be read
    [aOnset,Analog] = Widefield_LoadData(cFile,'Analog'); %load analog data
    triggerTime = aOnset(1) * (86400*1e3) - 50; Analog = double(Analog); Analog(2,:) = smooth(Analog(2,:));
    
    %% check for missing frames or dropped frames in the camera and throw warning if anything seems off
%     dSwitch = (diff((Analog(2,:)./std(Analog(2,:)) > 1)) == 1) + (diff((Analog(2,:)./std(Analog(2,:)) > 1)) == -1); %times when the photodiode detects a frame change
%     dSwitch(dSwitch<0) = 0;
%     if sum(diff(find(dSwitch)) > StimData.sRate*1000*1.5) > 1 %if time difference between two frames is larger as 1.5 times the display refresh rate
%         warning([num2str(sum(diff(find(dSwitch)) > StimData.sRate*1000*1.5)) ' skipped frames in trial ' int2str(iTrials) '; based on photodiode']);
%         warning(['Average time difference between stimulation frames: ' num2str(mean(diff(find(dSwitch)))) ' ms']);
%     end
%     firstFrame = triggerTime(1) + find(dSwitch,1); %timePoint of stimulus onset, based on photodiode signal.

    dSwitch = (diff(Analog(3,:) > 1000) == 1) + (diff(Analog(3,:) < 1000) == 1); %times when the stimulator triggers a frame change
    dSwitch(dSwitch<0) = 0;
    if any(diff(find(dSwitch)) > StimData.sRate*1000*1.5) %if time difference between two frames is larger as 1.5 times the display refresh rate
        warning([num2str(sum(diff(find(dSwitch)) > StimData.sRate*1000*1.5)) ' skipped frames in trial ' int2str(iTrials) '; based on trigger signal']);
        warning(['Average time difference between stimulation frames: ' num2str(mean(diff(find(dSwitch)))) ' ms']);
    end
    firstFrame = triggerTime(1) + find(dSwitch,1); %timePoint of stimulus onset, based on trigger signal.

    sFrameDur = round(mean(diff(frameTimes - frameTimes(1)))); %average duration between acquired frames
    dSwitch = diff(frameTimes - frameTimes(1)); %duration between all acquired frames
    if any(abs(dSwitch - sFrameDur) > sFrameDur*1.5)
        warning([num2str(sum(abs(dSwitch - sFrameDur) > sFrameDur*1.5)) ' lost camera frames in trial ' int2str(iTrials) '; something broken with camera settings?']);
        warning(['Average time difference between acquired frames: ' num2str(sFrameDur) ' ms']);
    end
    
    %% compute the amount of required frames and collect from data
    stimFrames = round(StimDur(iTrials)*1000 / sFrameDur); %required frames for stim sequence
    ind = find(diff((frameTimes-firstFrame) <= 0)); %index of frame that is closest to stimulus onset
    numCycles = StimDur(iTrials)/(1/barFreq(iTrials)); %number of cycles in current trial
    cycleFrames = round(stimFrames/numCycles); %frames per cycle

    Data(:,:,1:ind) = []; %throw away baseline
    Data = Data(:,:,1:stimFrames); %only use stimulation part of the data;

    binSize = floor((max(size(Data(:,:,1)))/winSize)/40); %compute binsize to get closest to 40 pixels/mm (recommended for segmentation code).
    if winSize > 1 && winSize < inf
        bData = arrayResize(Data,binSize); %do spatial binning
    else
        bData = Data;
    end
    clear Data

    %% running average
    for x = 1:length(nTrials)
        cTrialCnt = rem(condCnt(iBarConds(:,iTrials))+1,nTrials(x)); %current trial cycle for running average (reset to 1, when 'nTrials' is reached)
        
        if cTrialCnt == 1 || condCnt(iBarConds(:,iTrials)) == 0 %start cycle for running average
            avgData{iBarConds(:,iTrials),x} = bData; %starting dataset for running average with set trialcount
        else
            avgData{iBarConds(:,iTrials),x} = (avgData{iBarConds(:,iTrials),x}.*cTrialCnt + bData) ./ (cTrialCnt+1); %produce running average
        end
        
        if cTrialCnt == 0 %reached requested trialcount. Compute fourier transform and increase counter
            temp = fft(avgData{iBarConds(:,iTrials),x},[],3);
            fTransform{iBarConds(:,iTrials),x}(TrialCnts(iBarConds(:,iTrials),x),:,:) = squeeze(temp(:,:,numCycles+1)); clear temp
            TrialCnts(iBarConds(:,iTrials),x) = TrialCnts(iBarConds(:,iTrials),x)+1; %increase counter for running average when required trialcount is reached
        end
    end
    condCnt(iBarConds(:,iTrials)) =  condCnt(iBarConds(:,iTrials))+1;
    disp(['Done loading trial ' int2str(iTrials) '/' int2str(max(Trials))]);
end
clear Data bData

%% do fft analysis to get phase and magnitude maps
StimDur = unique(StimDur); %Duration of a given trial
barFreq = unique(barFreq); %Duration of a given trial
screenSize = fliplr(textscan(StimData.handles.ScreenSizeAngle,'%f%c%f')); %get screen size in visual angles. Flip this so first entry is for horizontal screen size, and third entry for vertical screen size.

% for iTrials = 1:length(nTrials)
for iTrials = 1
    Cnt = 1;
    for iConds = [1 3]; %this expects 4 directions to construct horizontal and vertical map
        for iRuns = 1:size(fTransform{iConds,iTrials},1)
        
           magMaps{Cnt,iRuns} = imrotate(squeeze(abs(fTransform{iConds,iTrials}(iRuns,:,:).*fTransform{iConds+1,iTrials}(iRuns,:,:))),rotateImage); %combined magnitude map.
           phaseMaps{Cnt,iRuns} = imrotate(squeeze(angle(fTransform{iConds,iTrials}(iRuns,:,:)./fTransform{iConds+1,iTrials}(iRuns,:,:))),rotateImage); %combined phase map (horizontal is iConds=1, vertical is iConds =3)
           phaseMaps{Cnt,iRuns}(isnan(phaseMaps{Cnt,iRuns}(:))) = 0;
           phaseMaps{Cnt,iRuns} = spatialFilterGaussian(phaseMaps{Cnt,iRuns}/pi*screenSize{iConds}/2,phaseMapSmth); % Translate to visual angles. Half of screenSize is used for each direction (assuming that animals eye is centered on the screen).

        end
        cPhaseMaps{Cnt,iTrials} = median(cat(3,phaseMaps{Cnt,:}),3);
        Cnt = Cnt+1;
    end
    cMagMaps{iTrials} = median(cat(3,magMaps{:}),3);
    cMagMaps{iTrials} =(cMagMaps{iTrials}-min(cMagMaps{iTrials}(:)))./(max(cMagMaps{iTrials}(:))- min(cMagMaps{iTrials}(:))); %normalize between 0 and 1

    % compute visual field sign maps. First compute gradients and atan - same as in the 2014 Callaway paper.
    for iRuns = 1:size(fTransform{iConds,iTrials},1)
        [dhdx, dhdy] = gradient(phaseMaps{1,iRuns});
        [dvdx, dvdy] = gradient(phaseMaps{2,iRuns});
        
        graddir_hor = atan2(dhdy,dhdx);
        graddir_vert = atan2(dvdy,dvdx);
        vdiff = exp(1i*graddir_hor) .* exp(-1i*graddir_vert);
        
        VFS{iRuns} = sin(angle(vdiff)); %Visual field sign map
        VFS{iRuns} = spatialFilterGaussian(VFS{iRuns},smth);
    end
    
    cVFS{1,iTrials} = median(cat(3,VFS{:}),3);
    clear magMaps phaseMaps VFS
    
    h = figure;
    subplot(2,2,1);
    imagesc(cPhaseMaps{1,iTrials});axis image; colormap hsv; colorbar; freezeColors;
    title(['Horizontal - nTrials = ' num2str(nTrials(iTrials))]);
    subplot(2,2,2);
    imagesc(cPhaseMaps{2,iTrials});axis image; colormap hsv; colorbar; freezeColors;
    title(['Vertical - nTrials = ' num2str(nTrials(iTrials))]);
    subplot(2,2,3);
    imagesc(cMagMaps{iTrials});axis image; colorbar; colormap jet;
    title('Mean Magnitude');
    subplot(2,2,4);
    imagesc(spatialFilterGaussian(cVFS{1,iTrials},smth)); axis image;colorbar
    title(['VisualFieldSign - binSize = ' num2str(binSize) '; smth = ' num2str(smth)]);
    savefig(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_allPlots_ ' int2str(nTrials(iTrials)) '_trials.fig']);
    h.PaperUnits = 'inches';
    set(h, 'PaperPosition', [0 0 15 15]);
    saveas(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_allPlots_ ' int2str(nTrials(iTrials)) '_trials.jpg'])
    clear h
end

trialSelect = 1; %use this to select which trialcount should be used for subsequent figures

%% create phase and magnitude maps for all requested trialcounts and show in seperate figures run visual segmentation code. 
% The last 2 inputs are the vfs threshold and the smoothing factor.  
% Seems like those needs to played with to produce a reasonable result.
% try
%     [im,nf_im] = getMouseAreasX(cPhaseMaps{1,trialSelect},cPhaseMaps{2,trialSelect},cMagMaps{trialSelect},40,4,3);
%     im = imresize(nf_im,binSize);  %use non-fused image
% catch ME
%     disp(['Error in ' ME.stack(1).name ': ' ME.message])
%     im = [];
% end
% im = imresize(im,binSize);
im = [];

%% get phase and amplitude + vessel map for plotting
plotPhaseMap = spatialFilterGaussian(cVFS{1,trialSelect},smth);
plotPhaseMap = imresize(plotPhaseMap,binSize);
plotPhaseMap =(plotPhaseMap-min(plotPhaseMap(:)))./(max(plotPhaseMap(:))- min(plotPhaseMap(:))); %normalize between 0 and 1

plotAmpMap = spatialFilterGaussian(imresize(cMagMaps{trialSelect},binSize),25); %smoothed magnitude map
plotAmpMap =(plotAmpMap-min(plotAmpMap(:)))./(max(plotAmpMap(:))- min(plotAmpMap(:))); %normalize between 0 and 1
plotAmpMap = spatialFilterGaussian(plotAmpMap,25); %smoothed magnitude map
% plotAmpMap =(plotAmpMap-min(plotAmpMap(:)))./(max(plotAmpMap(:))- min(plotAmpMap(:))); %normalize between 0 and 1

% cFile = [path Animal '\PhaseMap\' Session '\Frames_' num2str(Trials(1)) '.dat']; %current file to be read
% [~,Data] = Widefield_LoadData(cFile,'Frames'); %load video data      
% snap = smooth2a(mean(Data,4),2,2); %use this to use raw data instead of green light vessel image
% snap = double(imrotate(snap,rotateImage));
% snap =(snap-min(snap(:)))./(max(snap(:))- min(snap(:))); %normalize between 0 and 1

cFile = ls([path Animal '\PhaseMap\' Session '\Snapshot_1.mat']); %get vessel image.
load([path Animal '\PhaseMap\' Session '\' cFile]);
snap = double(imrotate(snap,rotateImage));
snap =(snap-min(snap(:)))./(max(snap(:))- min(snap(:))); %normalize between 0 and 1

if size(snap,1) > size(plotPhaseMap,1) %resize if vessel map is larger
    plotPhaseMap = imresize(plotPhaseMap,size(snap,1)/size(plotPhaseMap,1));
    plotAmpMap = imresize(plotAmpMap,size(snap,1)/size(plotAmpMap,1));
end

%% plot vessel map with overlayed sign map
h = figure;
imagesc(plotPhaseMap); colormap jet; 
caxis([0 1]);
title([Animal ' - PhaseMap - Color']);axis image
savefig(h,[path Animal '\PhaseMap\' Session '\' Animal '_RawPhaseMap.fig']);
saveas(h,[path Animal '\PhaseMap\' Session '\' Animal '_RawPhaseMap.jpg']);
save([path Animal '\PhaseMap\' Session '\plotPhaseMap.mat'],'plotPhaseMap');
save([path Animal '\PhaseMap\' Session '\plotAmpMap.mat'],'plotAmpMap');
save([path Animal '\PhaseMap\' Session '\segmentImg.mat'],'im');
save([path Animal '\PhaseMap\' Session '\cMagMaps.mat'],'cMagMaps');
save([path Animal '\PhaseMap\' Session '\cPhaseMaps.mat'],'cPhaseMaps');
save([path Animal '\PhaseMap\' Session '\cVFS.mat'],'cVFS');
close(h);

h = figure;
imagesc(snap);axis image; colormap gray; freezeColors; hold on
vfsIm = imagesc(plotPhaseMap); colormap jet; 
caxis([0 1]);
set(vfsIm,'AlphaData',plotAmpMap); axis image
title([Animal ' - PhaseMap - Vesselmap + Color'])
savefig(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap.fig']);
saveas(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap.jpg'])
close(h);

h = figure;
imagesc(snap);axis image; colormap gray; freezeColors; hold on
contour(im,[.5 .5],'w','linewidth',3); axis image
title([Animal ' - PhaseMap - Vesselmap + Outline'])
savefig(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_outlined_noColor.fig']);
saveas(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_outlined_noColor.jpg'])
close(h);

h = figure;
imagesc(snap);axis image; colormap gray; freezeColors; hold on
vfsIm = imagesc(plotPhaseMap); colormap jet; 
caxis([0 1]);
set(vfsIm,'AlphaData',plotAmpMap); axis image
contour(im,[.5 .5],'w','linewidth',3); axis image
title([Animal ' - PhaseMap - Vesselmap + Color + Outline'])
savefig(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_outlined.fig']);
saveas(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_outlined.jpg'])


%% show stereotactic coordinates
if ~isempty(pxPerMM) %this can be used to compute stereotactic coordinates when supplying the resolution in pixels/mm.
    h1 = figure;
    Check = false;
    cFile = [path Animal '\PhaseMap\' Session '\' fName num2str(Trials(1)) '.dat']; %current file to be read
    [~,Data] = Widefield_LoadData(cFile,'Frames'); %load video data
        
    while ~Check
        
        cFile = ls([path Animal '\PhaseMap\' Session '\Snapshot_1.mat']); %get vessel image.
        load([path Animal '\PhaseMap\' Session '\' cFile]);
%         snap = mean(Data,4); %use this to use raw data instead of green light vessel image
        snap = double(imrotate(snap,rotateImage)); %this should be the green light vessel image
        snap =(snap-min(snap(:)))./(max(snap(:))- min(snap(:))); %normalize between 0 and 1
        
        plotAmpMap = imresize(cMagMaps{trialSelect},binSize); %smoothed magnitude map
        plotAmpMap = spatialFilterGaussian(plotAmpMap,50); %smoothed magnitude map
        plotAmpMap =(plotAmpMap-min(plotAmpMap(:)))./(max(plotAmpMap(:))- min(plotAmpMap(:))); %normalize between 0 and 1
        plotAmpMap(:,end/2:end) =((temp-min(temp(:)))./(max(temp(:))- min(temp(:))))/2; %normalize between 0 and 1
        plotAmpMap = spatialFilterGaussian(plotAmpMap,25); %smoothed magnitude map
        
        plotPhaseMap = spatialFilterGaussian(cVFS{1,trialSelect},smth);
        plotPhaseMap = imresize(plotPhaseMap,binSize);
        plotPhaseMap =(plotPhaseMap-min(plotPhaseMap(:)))./(max(plotPhaseMap(:))- min(plotPhaseMap(:))); %normalize between 0 and 1
        
        if size(snap,1) > size(plotPhaseMap,1) %resize if vessel map is larger
            plotPhaseMap = imresize(plotPhaseMap,size(snap,1)/size(plotPhaseMap,1));
            plotAmpMap = imresize(plotAmpMap,size(snap,1)/size(plotAmpMap,1));
        end
        
        imagesc(snap);axis image; colormap gray; caxis([0 1]);set(gca,'linewidth',1)
        grid(gca,'on');grid minor;set(gca,'GridColor','w');
        set(gca,'xTick',1:pxPerMM:size(snap,1))
        set(gca,'yTick',1:pxPerMM:size(snap,2))
        
        Wait = input('Roate vesselpic(deg).\n ','S');
        if ~isnan(str2double(Wait))
            fineRotate = str2double(Wait);
            snap = imrotate(snap,fineRotate);
            plotPhaseMap = imrotate(plotPhaseMap,fineRotate);
            plotAmpMap = imrotate(plotAmpMap,fineRotate);
            
            imagesc(snap);axis image; colormap gray; caxis([0 1]);
            grid(gca,'on');grid minor;set(gca,'GridColor','w');
            set(gca,'xTick',1:pxPerMM:size(snap,1))
            set(gca,'yTick',1:pxPerMM:size(snap,2))
        end
        
        Wait = input('Type "y" to continue or any other key to set angle again. \n ','S');
        if strcmpi(Wait,'y')
            Check = true;grid minor;
            set(gca,'linewidth',3)
            hold on
            freezeColors;
        end
    end
    
    %shift grid in x
    Check = false;
    while ~Check
        title('Select bregma coordinates to align grid');
        [x,y] = ginput(1);
        xVec = x-floor(x/pxPerMM)*pxPerMM:pxPerMM:size(snap,1);
        xLabel = num2str(abs((1:length(xVec))-ceil(x/pxPerMM))');
        set(gca,'xTick',xVec); set(gca,'xTickLabel',xLabel)
        
        yVec = y-floor(y/pxPerMM)*pxPerMM:pxPerMM:size(snap,2);
        yLabel = num2str(abs((1:length(yVec))-ceil(y/pxPerMM))');
        set(gca,'yTick',yVec); set(gca,'yTickLabel',yLabel)
        
        Wait = input('Type "y" to continue or any other key to set bregma again. \n ','S');
        if strcmpi(Wait,'y')
            Check = true;
        end
    end
    title('Vessel image; FFT projection + stereotactic coordinates');
    xlabel('Mediolateral from Bregma(mm)');
    ylabel('Anterioposterior from Bregma (mm)');
    vfsIm = imagesc(plotPhaseMap); colormap jet; caxis([0.25 0.75]);
    set(vfsIm,'AlphaData',plotAmpMap); axis image
    
    h = figure;
    subplot(1,2,1)
    imagesc(snap);axis image; colormap gray; caxis([0 1]);
    grid(gca,'on');set(gca,'GridColor','w');freezeColors
    axis image
    title('Vessel image');
    set(gca,'xTick',1:pxPerMM:size(snap,1));
    set(gca,'xTick',xVec); set(gca,'xTickLabel',xLabel)
    set(gca,'yTick',1:pxPerMM:size(snap,2));        
    set(gca,'yTick',yVec); set(gca,'yTickLabel',yLabel)
    xlabel('Mediolateral from Bregma(mm)');
    ylabel('Anterioposterior from Bregma (mm)');

    subplot(1,2,2)
    imagesc(plotPhaseMap); colormap jet; 
    caxis([0.25 0.75]);
    grid(gca,'on');set(gca,'GridColor','k');
    axis image; colormap jet
    title('Visual field sign');
    set(gca,'xTick',1:pxPerMM:size(snap,1));
    set(gca,'xTick',xVec); set(gca,'xTickLabel',xLabel)
    set(gca,'yTick',1:pxPerMM:size(snap,2));        
    set(gca,'yTick',yVec); set(gca,'yTickLabel',yLabel)
    xlabel('Mediolateral from Bregma(mm)');
    ylabel('Anterioposterior from Bregma (mm)');
    h.PaperUnits = 'inches';
    set(h, 'PaperPosition', [0 0 15 15]);

end

%embedded some functions to make the code work as stand-alone
function img = spatialFilterGaussian(img, sigma)
if sigma > 0 && (numel(img) ~=  sum(sum(isnan(img))))
    hh = fspecial('gaussian',size(img),sigma);
    hh = hh/sum(hh(:));
    img = ifft2(fft2(img).*abs(fft2(hh)));
end

function [header,data] = Widefield_LoadData(cPath,Condition,varargin)
% short routine to load data from WidefieldImager code.
% cPath is the path of the file that should be opened. Condition is the
% type of data file which will determine the way the data file is read.
% Optional input 'pInd' defines a single pixel from which a data vector should
% be extracted. pInd is a two-value vector for x-y pixel coordinates. This
% means X and Y for image coordinates NOT matlab coordinates (which are
% usually inverse).

if ~isempty(varargin)
    pInd = varargin{1}; %index for pixel to be extracted
    if length(pInd) ~= 2 || ~isnumeric(pInd)
        error('Invalid input for index of selected pixel')
    end
end

fID = fopen(cPath);
switch lower(Condition)
    case 'analog'
        hSize = fread(fID,1,'double'); %header size
        header = fread(fID,hSize,'double'); %Metadata. Default is: 1 = Time of Acquisition onset, 2 = Number of channels, 3 = number of values per channel
        data = fread(fID,[header(end-1),header(end)],'uint16=>uint16'); %get data. Last 2 header values should contain the size of the data array.
    case 'frames'
        hSize = fread(fID,1,'double'); %header size
        header = fread(fID,hSize,'double'); %Metadata. Default is: 1:x = Absolute timestamps for each frame, Last 4 values: Size of each dimensions in the matrix
        
        if ~isempty(varargin) %if extracting single pixel information
            imSize = (header(find(diff(header) < -1e4) + 1)*header(find(diff(header) < -1e4) + 2))-1; %number of datapoints to make a single image minus one. skip that many datapoints to stick to the same pixel when using fread.
            imStart = ((pInd(1)-1)*header(find(diff(header) < -1e4) + 1))+pInd(2)-1; %first value for selected pixel
            fseek(fID,imStart*2,'cof'); %shift file pointer to the right pixel to start data extraction from file
            data = fread(fID,header(find(diff(header) < -1e4) + 4),'uint16=>uint16',imSize*2); %get data.
            if length(data) ~= header(end)
                error('Could not extract all data values from pixel')
            end
        else
            data = fread(fID,[prod(header(find(diff(header) < -1e4) + 1 : end)),1],'uint16=>uint16'); %get data. Last 4 header values should contain the size of the data array.
            if length(data) ~= prod(header(find(diff(header) < -1e4) + 1 : end)) %if insufficient data is found in .dat file. Sometimes fread does not get all values from file when reading from server.
                fclose(fID);fID = fopen(cPath); %try loading data again
                hSize = fread(fID,1,'double'); %header size
                header = fread(fID,hSize,'double'); %Metadata. Defautlt is: 1:x = Absolute timestamps for each frame, Last 4 values: Size of each dimensions in the matrix
                data = fread(fID,[prod(header(find(diff(header) < -1e4) + 1 : end)),1],'uint16=>uint16'); %get data. Last 4 header values should contain the size of the data array.
            end
            data = reshape(data,header(find(diff(header) < -1e4) + 1 : end)'); %reshape data into matrix
        end
end
fclose(fID);

function frameRecon = svdFrameReconstruct(U, V)
% function frameRecon = svdFrameReconstruct(U, V)
% U is Y x X x nSVD
% V is nSVD x nFrames

% reshape U to be nPix x nSVD
Ur = reshape(U, size(U,1)*size(U,2), size(U,3));

% multiply and reshape back into Y x X
frameRecon = reshape(Ur*V, size(U,1), size(U,2), size(V,2));