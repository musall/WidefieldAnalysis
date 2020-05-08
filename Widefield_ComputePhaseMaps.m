function Widefield_ComputePhaseMaps(Animal,Session,winSize,smth,rotateImage,nTrials,pxPerMM)
dbstop Widefield_ComputePhaseMaps 221;

%% Set basic variables
fclose('all');
path = '\\grid-hs\churchland_nlsas_data\data\WidefieldImager\Animals\';   %Widefield data path
% path = '\\grid-hs\churchland_hpc_home\smusall\WidefieldImager\Animals\';   %Widefield data path
% path = 'W:\data\WidefieldImager\Animals\';   %Widefield data path
% path = 'X:\smusall\WidefieldImager\Animals\';   %Widefield data path
% path = 'Q:\';   %Widefield data path
checkWheel = false;                      %check for wheel motion. This assumes velocity data in the second last analog channel and will try to find a running area.
phaseMapSmth = 2;                       % smoothing factor for phaseMap
opts.fName = 'Frames';
opts.fPath = [path Animal filesep 'PhaseMap' filesep Session filesep];
opts.plotChans = true;
opts.trigLine = [6 7];
% opts.trigLine = NaN;
opts.blueHigh = true;
opts.stimLine = 3;
opts.preStim = 1;
opts.alignRes = 10;
opts.hemoCorrect = true;
opts.frameRate = 30; %framerate in Hz
opts.fileExt = 'mp4'; %type of video file. Default is mj2.

if ~exist('winSize','var')
    winSize = 10;                        %size of the field in mm. This is to determine the spatial binning to get as close to 40pix/mm as possible. This is advised for the visual segmentation code.
end

if ~exist('smth','var')
    smth = 2;                         %smoothing factor for gaussian smooth on maps
end

if ~exist('rotateImage','var')
    rotateImage = 0;                    %rotate image if the orientation is not as required
end

if ~exist('nTrials','var')
    nTrials = 1;                        %number of trials for phasemaps. If nCycles <=0 or > the number of cyles in the data. It will use all available cycles for a single phasemap.
end

if ~exist('pxPerMM','var')
    pxPerMM = [];                       %mapping from pixels to real space in pixels/mm. Typically, this is 165px/mm if fully zoomed in and 128px/mm if fully zoomed out. See the RulerPics in the WidefieldImager folder to check for yourself. 
end


%% data source
disp(['Current path: ' path Animal '\PhaseMap\' Session]); tic

%get overview of trials by looking at analog data
useVideo = false;
temp = ls([opts.fPath opts.fName '*dat']);
temp1 = reshape(temp',1,numel(temp));
temp1 = strrep(temp1,'.dat','    ');
if isempty(temp1)
    temp = ls([opts.fPath opts.fName '*' opts.fileExt]); 
    temp1 = reshape(temp',1,numel(temp));
    temp1 = strrep(temp1,['.' opts.fileExt],'    ');
    useVideo = true;
end
temp = reshape(temp1,size(temp'))';clear temp1
Trials = sort(str2num(temp(:,length(opts.fName)+2:end)));clear temp %identified trials

if rem(length(Trials),4) > 0 %trials has to be a divider of 4
    Trials(end-rem(length(Trials),4)+1) = [];
end

cFile = ls([opts.fPath Animal '*settings.mat']);
load([path Animal '\PhaseMap\' Session '\' cFile]);

if length(Trials) ~= str2num(StimData.handles.NrTrials) %compare recorded data to set amount of trials
    warning(['Recorded trials (' num2str(length(Trials)) ') are unequal to settings in visual stimulator (' StimData.handles.NrTrials ')'])
end

%get different bar direction and orentiations for individual trials
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
opts.postStim = StimDur(1);

% check for number of cycles and trials in dataset
if length(numCycles) ~= 1
    error('Number of cycles per trial is inconsistent. This code is not meant to handle that.')
end
nTrials(nTrials <= 0 | nTrials > length(Trials)/4) = length(Trials)/4; 
nTrials = [nTrials floor(length(Trials)/4)]; nTrials = unique(nTrials);
disp(['Trials per condition: ' num2str(length(Trials)/4) ' - Computing phasemaps from [' num2str(nTrials) '] trials']);
TrialCnts = ones(4,length(nTrials));
avgData = cell(4,length(nTrials)); % averaged sequence for each condition and trialcount
fTransform = cell(4,length(nTrials)); % fourier transforms for each condition and trialcount 

%% get single trials and average for each condition to end up with one sequence
tic;
condCnt = zeros(1,4); %counter for how many sequences were saved in each condition
for iTrials = Trials'
    %% load data
    %     Data = Widefield_HemoCorrect(blueData,hemoData,1:ceil(opts.preStim * sRate),5); %hemodynamic correction
    
    if strcmpi(opts.fName,'acFrames')
        cFile = [opts.fPath filesep opts.fName '_' num2str(iTrials) '.dat']; %current file to be read
        [header,Data] = Widefield_LoadData(cFile,'Frames'); %load video data
        frameTimes = header(1:end-ndims(Data)); %extract frame times from header and convert to millisecond timestamps
        Data(:,:,1:opts.preStim * round(opts.frameRate) - 1) = []; %throw away baseline
            
        binSize = floor((max(size(Data(:,:,1)))/winSize)/40); %compute binsize to get closest to 40 pixels/mm (recommended for segmentation code).
        if winSize > 1 && winSize < inf
            Data = arrayResize(Data,binSize); %do spatial binning
        end
    else
        if useVideo
            [blueData,frameTimes,hemoData,~,~,~,~] = Widefield_SplitVideo(opts,iTrials,opts.fileExt); %load mixed or blue only data from video files
        else
            [blueData,frameTimes,hemoData,~,~,~,~] = Widefield_SplitChannels(opts,iTrials); %load mixed or blue only data from .dat files
        end
        toc;
        binSize = floor((max(size(blueData(:,:,1)))/winSize)/40); %compute binsize to get closest to 40 pixels/mm (recommended for segmentation code).
        if winSize > 1 && winSize < inf
            blueData = arrayResize(blueData,binSize); %do spatial binning
            hemoData = arrayResize(hemoData,binSize); %do spatial binning
        end
        if opts.hemoCorrect
            Data = Widefield_HemoCorrect(blueData,hemoData,1:ceil(opts.preStim * opts.frameRate),5); %hemodynamic correction
        else
            Data = blueData; clear blueData
        end
        Data(:,:,1:ceil(opts.preStim * opts.frameRate)) = []; %throw away baseline
        toc;
    end
%     meanTrace(:,iTrials) = mean(reshape(Data,[],size(Data,3))); %get a mean average for the whole frame to check for suspicous activity

    cFile = [opts.fPath filesep 'Analog_' num2str(iTrials) '.dat']; %current file to be read
    [~,Analog] = Widefield_LoadData(cFile,'Analog'); %load analog data
    Analog = double(Analog); Analog(2,:) = smooth(Analog(2,:));
    
    %% check for missing frames or dropped frames in the camera and throw warning if anything seems off
    dSwitch = (diff((Analog(2,:)./std(Analog(2,:)) > 1)) == 1) + (diff((Analog(2,:)./std(Analog(2,:)) > 1)) == -1); %times when the photodiode detects a frame change
    dSwitch(dSwitch<0) = 0;
    if sum(diff(find(dSwitch)) > StimData.sRate*1000*1.5) > 1 %if time difference between two frames is larger as 1.5 times the display refresh rate
        warning([num2str(sum(diff(find(dSwitch)) > StimData.sRate*1000*1.5)) ' skipped frames in trial ' int2str(iTrials) '; based on photodiode']);
        warning(['Average time difference between stimulation frames: ' num2str(mean(diff(find(dSwitch)))) ' ms']);
    end

    dSwitch = (diff(Analog(3,:) > 1000) == 1) + (diff(Analog(3,:) < 1000) == 1); %times when the stimulator triggers a frame change
    dSwitch(dSwitch<0) = 0;
    if any(diff(find(dSwitch)) > StimData.sRate*1000*1.5) %if time difference between two frames is larger as 1.5 times the display refresh rate
        disp(['Warning: ' num2str(sum(diff(find(dSwitch)) > StimData.sRate*1000*1.5)) ' skipped frames in trial ' int2str(iTrials) '; based on trigger signal']);
        disp(['Average time difference between stimulation frames: ' num2str(mean(diff(find(dSwitch)))) ' ms']);
    end

    sFrameDur = round(1000/opts.frameRate); %average duration between acquired frames
%     sFrameDur = round(mean(diff(frameTimes - frameTimes(1)))); %average duration between acquired frames
    dSwitch = diff(frameTimes - frameTimes(1)); %duration between all acquired frames
    if any(abs(dSwitch - sFrameDur) > sFrameDur*1.5)
        warning([num2str(sum(abs(dSwitch - sFrameDur) > sFrameDur*1.5)) ' lost camera frames in trial ' int2str(iTrials) '; something broken with camera settings?']);
        warning(['Average time difference between acquired frames: ' num2str(sFrameDur) ' ms']);
    end
    
    %% compute the amount of required frames and collect from data
    if numCycles ~= floor(((1/opts.frameRate)*size(Data,3))/(1/barFreq(iTrials))) %make sure that number of cycles matches the current data
        error('Number of presented cycles does not match the length of available data')
    end
    
    %% running average
    for x = 1:length(nTrials)
        cTrialCnt = rem(condCnt(iBarConds(:,iTrials))+1,nTrials(x)); %current trial cycle for running average (resets to 1, after 'nTrials' is reached)
               
        if (cTrialCnt == 1 && condCnt(iBarConds(:,iTrials)) >= nTrials(x)) || condCnt(iBarConds(:,iTrials)) == 0 || nTrials(x) == 1 %start cycle for running average
            avgData{iBarConds(:,iTrials),x} = Data; %starting dataset for running average with set trialcount
        else
            avgData{iBarConds(:,iTrials),x} = (avgData{iBarConds(:,iTrials),x}.*condCnt(iBarConds(:,iTrials)) + Data) ./ (condCnt(iBarConds(:,iTrials))+1); %produce running average
        end
        
        if cTrialCnt == 0 %reached requested trialcount. Compute fourier transform and increase counter
            temp = fft(avgData{iBarConds(:,iTrials),x},[],3);
            fTransform{iBarConds(:,iTrials),x}(TrialCnts(iBarConds(:,iTrials),x),:,:) = squeeze(temp(:,:,numCycles + 1)); clear temp
            TrialCnts(iBarConds(:,iTrials),x) = TrialCnts(iBarConds(:,iTrials),x)+1; %increase counter for running average when required trialcount is reached
        end
    end
    condCnt(iBarConds(:,iTrials)) =  condCnt(iBarConds(:,iTrials))+1;
    disp(['Done loading trial ' int2str(iTrials) '/' int2str(max(Trials))]);
end
clear Data

%% plot activity traces across trials
% figure
% plot(meanTrace,'linewidth',1); axis square;
% xlabel('Frames'); ylabel('Mean signal change');

%% get vessel image.
cFile = ls([opts.fPath 'Snapshot_1.mat']); 
load([opts.fPath cFile], 'snap');
snap = double(imrotate(snap,rotateImage));
snap =(snap-min(snap(:)))./(max(snap(:))- min(snap(:))); %normalize between 0 and 1
snap = snap(1:size(avgData{1},1),1:size(avgData{1},2));

% if winSize > 1 && winSize < inf
%     snap = arrayResize(snap,binSize); %do spatial binning
% end
% snap(snap>.25) = .25;

% cFile = [path Animal '\PhaseMap\' Session '\Frames_' num2str(Trials(1)) '.dat']; %current file to be read
% [~,Data] = Widefield_LoadData(cFile,'Frames'); %load video data      
% snap = smooth2a(mean(Data,4),2,2); %use this to use raw data instead of green light vessel image
% snap = double(imrotate(snap,rotateImage));
% snap =(snap-min(snap(:)))./(max(snap(:))- min(snap(:))); %normalize between 0 and 1


%% do fft analysis to get phase and magnitude maps
StimDur = unique(StimDur); %Duration of a given trial
barFreq = unique(barFreq); %Duration of a given trial
screenSize = fliplr(textscan(StimData.handles.ScreenSizeAngle,'%f%c%f')); %get screen size in visual angles. Flip this so first entry is for horizontal screen size, and third entry for vertical screen size.
nPairs = min(cellfun(@(x) size(x,1), fTransform(:,1)));
clear phaseMaps magMaps
trialSelect = 1; %use this to select which trialcount should be used for subsequent figures

% for iTrials = 1:length(nTrials)
for iTrials = trialSelect
    Cnt = 1;
    for iConds = [1 3] %this expects 4 directions to construct horizontal and vertical map
        for iRuns = 1:nPairs
            
           magMaps{Cnt,iRuns} = imrotate(squeeze(abs(fTransform{iConds,iTrials}(iRuns,:,:).*fTransform{iConds+1,iTrials}(iRuns,:,:))),rotateImage); %combined magnitude map.
           
           a1 = mod(-angle(fTransform{iConds,iTrials}(iRuns,:,:)), 2 * pi);
           a2 = mod(-angle(fTransform{iConds+1,iTrials}(iRuns,:,:)), 2 * pi);
           a1(isnan(a1)) = 0; a2(isnan(a2)) = 0;
           phaseMaps{Cnt,iRuns} = imrotate(squeeze((a2 - a1) / 2),rotateImage);
           phaseMaps{Cnt,iRuns} = spatialFilterGaussian(phaseMaps{Cnt,iRuns}*screenSize{iConds}/2,phaseMapSmth); % Translate to visual angles. Half of screenSize is used for each direction (assuming that animals eye is centered on the screen).

%            phaseMaps{Cnt,iRuns} = imrotate(squeeze(angle(fTransform{iConds,iTrials}(iRuns,:,:)./fTransform{iConds+1,iTrials}(iRuns,:,:))),rotateImage); %combined phase map (horizontal is iConds=1, vertical is iConds =3)
% %            phaseMaps{Cnt,iRuns} = mod(phaseMaps{Cnt,iRuns},2*pi); %this step is only required when not dividing in fourier space when computing phase maps. Remaps angles to -pi:pi instead of 0:-pi:pi:0 which breacks the sign map.
%            phaseMaps{Cnt,iRuns}(isnan(phaseMaps{Cnt,iRuns}(:))) = 0;
%            phaseMaps{Cnt,iRuns} = spatialFilterGaussian(phaseMaps{Cnt,iRuns}/pi*screenSize{iConds}/2,phaseMapSmth); % Translate to visual angles. Half of screenSize is used for each direction (assuming that animals eye is centered on the screen).

        end
        cPhaseMaps{Cnt,iTrials} = median(cat(3,phaseMaps{Cnt,:}),3);
        Cnt = Cnt+1;
    end
    cMagMaps{iTrials} = median(cat(3,magMaps{:}),3);
%     cMagMaps{iTrials}(snap < prctile(snap(:),75)) = 0;
    cMagMaps{iTrials} =(cMagMaps{iTrials}-min(cMagMaps{iTrials}(:)))./(max(cMagMaps{iTrials}(:))- min(cMagMaps{iTrials}(:))); %normalize between 0 and 1

    % compute visual field sign maps
    for iRuns = 1:nPairs
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
%     savefig(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_allPlots_ ' int2str(nTrials(iTrials)) '_trials.fig']);
%     h.PaperUnits = 'inches';
%     set(h, 'PaperPosition', [0 0 15 15]);
%     print(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_allPlots_ ' int2str(nTrials(iTrials)) '_trials.jpg'],'-djpeg')
%     clear h
end


%% get phase and amplitude + vessel map for plotting
plotPhaseMap = spatialFilterGaussian(cVFS{1,trialSelect},smth);
plotPhaseMap = imresize(plotPhaseMap,binSize);
plotPhaseMap =(plotPhaseMap-min(plotPhaseMap(:)))./(max(plotPhaseMap(:))- min(plotPhaseMap(:))); %normalize between 0 and 1

plotAmpMap = imresize(cMagMaps{trialSelect},binSize); %smoothed magnitude map
temp = plotAmpMap(:,1:end/2);
plotAmpMap(:,1:end/2) =(temp-min(temp(:)))./(max(temp(:))- min(temp(:))); %normalize between 0 and 1
temp = plotAmpMap(:,end/2:end);
plotAmpMap(:,end/2:end) =(temp-min(temp(:)))./(max(temp(:))- min(temp(:))); %normalize between 0 and 1
plotAmpMap = spatialFilterGaussian(plotAmpMap,25); %smoothed magnitude map
plotAmpMap =(plotAmpMap-min(plotAmpMap(:)))./(max(plotAmpMap(:))- min(plotAmpMap(:))); %normalize between 0 and 1

if size(snap,1) > size(plotPhaseMap,1) %resize if vessel map is larger
    plotPhaseMap = imresize(plotPhaseMap,size(snap,1)/size(plotPhaseMap,1));
    plotAmpMap = imresize(plotAmpMap,size(snap,1)/size(plotAmpMap,1));
end

%% create phase and magnitude maps for all requested trialcounts and show in seperate figures run visual segmentation code. 
% The last 2 inputs are the vfs threshold and the smoothing factor.  
% Seems like those needs to played with to produce a reasonable result.
try
    [im,nf_im] = getMouseAreasX(cPhaseMaps{1,trialSelect},cPhaseMaps{2,trialSelect},cMagMaps{trialSelect},40,4,3);
    im = imresize(nf_im,size(snap,1)/size(nf_im,1));  %use non-fused image
catch ME
    disp(['Error in ' ME.stack(1).name ': ' ME.message])
    im = [];
end
% im = imresize(im,binSize);

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
set(vfsIm,'AlphaData',plotAmpMap*1.5); axis image
title([Animal ' - PhaseMap - Vesselmap + Color'])
savefig(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap.fig']);
print(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap.jpg'],'-djpeg')
close(h);

% h = figure;
% imagesc(snap);axis image; colormap gray; freezeColors; hold on
% contour(im,[.5 .5],'w','linewidth',3); axis image
% title([Animal ' - PhaseMap - Vesselmap + Outline'])
% savefig(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_outlined_noColor.fig']);
% print(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_outlined_noColor.jpg'],'-djpeg')
% close(h);

h = figure;
imagesc(snap);axis image; colormap gray; freezeColors; hold on
vfsIm = imagesc(plotPhaseMap); colormap jet; 
caxis([0 1]);
set(vfsIm,'AlphaData',plotAmpMap); axis image
contour(im,[.5 .5],'w','linewidth',3); axis image
title([Animal ' - PhaseMap - Vesselmap + Color + Outline'])
savefig(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_outlined.fig']);
print(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_outlined.jpg'],'-djpeg')


%% show stereotactic coordinates
if ~isempty(pxPerMM)
    h1 = figure;
    Check = false;
%     cFile = [path Animal '\PhaseMap\' Session '\' fName num2str(Trials(1)) '.dat']; %current file to be read
%     [~,Data] = Widefield_LoadData(cFile,'Frames'); %load video data
        
    while ~Check
        
        cFile = ls([path Animal '\PhaseMap\' Session '\Snapshot_1.mat']); %get vessel image.
        load([path Animal '\PhaseMap\' Session '\' cFile]);
%         load('U:\space_managed_data\WidefieldImager\NewOptics\InfinityHigh_3ftLow.mat'); % ruler image
%         snap = mean(Data,4); %use this to use raw data instead of green light vessel image
        snap = double(imrotate(snap,rotateImage)); %this should be the green light vessel image
        snap =(snap-min(snap(:)))./(max(snap(:))- min(snap(:))); %normalize between 0 and 1
        
        plotAmpMap = imresize(cMagMaps{trialSelect},binSize); %smoothed magnitude map
        %         plotAmpMap = spatialFilterGaussian(plotAmpMap,50); %smoothed magnitude map
        %         plotAmpMap =(plotAmpMap-min(plotAmpMap(:)))./(max(plotAmpMap(:))- min(plotAmpMap(:))); %normalize between 0 and 1
        temp = plotAmpMap(:,1:end/2);
        plotAmpMap(:,1:end/2) =((temp-min(temp(:)))./(max(temp(:))- min(temp(:))))*1.25; %normalize between 0 and 1
        temp = plotAmpMap(:,end/2:end);
        plotAmpMap(:,end/2:end) =((temp-min(temp(:)))./(max(temp(:))- min(temp(:))))/2; %normalize between 0 and 1
        plotAmpMap = spatialFilterGaussian(plotAmpMap,25); %smoothed magnitude map
        
        plotPhaseMap = spatialFilterGaussian(cVFS{1,trialSelect},smth);
        plotPhaseMap = -imresize(plotPhaseMap,binSize);
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
        xVec = x-floor(x/pxPerMM)*pxPerMM:pxPerMM:size(snap,2);
        xLabel = num2str(abs((1:length(xVec))-ceil(x/pxPerMM))');
        set(gca,'xTick',xVec); set(gca,'xTickLabel',xLabel)
        
        yVec = y-floor(y/pxPerMM)*pxPerMM:pxPerMM:size(snap,1);
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
    vfsIm = imagesc(plotPhaseMap); colormap jet; caxis([0 1]);
    set(vfsIm,'AlphaData',plotAmpMap); axis image
    if ~checkWheel
        savefig(gcf,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_measureGridline.fig']);
        saveas(gcf,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_measureGridline.jpg'])
    end
    
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
    caxis([0 1]);
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
    if ~checkWheel
        savefig(gcf,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_GridlineCompare.fig']);
        saveas(gcf,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_GridlineCompare.jpg'])
    end
end

%% get wheel motion data
if checkWheel
    [aInd,runImg] = WidefieldImager_RunningMouse_v1(Animal,Session,binSize,smth,80,rotateImage+fineRotate);
    runImg(runImg==0) = NaN;
    runImg = imresize(runImg,binSize);
    aInd = aInd*binSize;
    
    %% plot motion data
    figure(h)
    subplot(1,3,3)
    imagesc(runImg); colormap jet;
    grid(gca,'on');set(gca,'GridColor','k');
    axis image; colormap jet
    title('Running mouse (deltaR/R)');
    set(gca,'xTick',1:pxPerMM:size(snap,1));
    set(gca,'xTick',xVec); set(gca,'xTickLabel',xLabel)
    set(gca,'yTick',1:pxPerMM:size(snap,2));        
    set(gca,'yTick',yVec); set(gca,'yTickLabel',yLabel)
    xlabel('Mediolateral from Bregma(mm)');
    ylabel('Anterioposterior from Bregma (mm)');hold on
    plot(aInd(1,:),aInd(2,:),'k','linewidth',2); axis image;
    savefig(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_GridlineCompare.fig']);
    saveas(h,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_GridlineCompare.jpg'])
    
    figure(h1)
    plot(aInd(1,:),aInd(2,:),'w','linewidth',2); axis image;
    savefig(h1,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_measureGridline.fig']);
    saveas(h1,[path Animal '\PhaseMap\' Session '\' Animal '_phaseMap_measureGridline.jpg'])
end

function img = spatialFilterGaussian(img, sigma)
if sigma > 0 && (numel(img) ~=  sum(sum(isnan(img))))
    hh = fspecial('gaussian',size(img),sigma);
    hh = hh/sum(hh(:));
    img = ifft2(fft2(img).*abs(fft2(hh)));
end
