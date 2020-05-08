load('C:\data\WidefieldImager\Animals\mSM21\MixedStim\17-Nov-2016\fftData.mat')
load('C:\data\WidefieldImager\Animals\mSM21\MixedStim\17-Nov-2016\avgTrace.mat')

%% normalize
for iCond = 1:length(avgTrace)
    ind = size(avgTrace{iCond},3);
    trace = squeeze(mean(avgTrace{iCond}(:,:,ind/4:ind/2),3));
    for iFrames = 1:size(avgTrace{iCond},3)
        avgTrace{iCond}(:,:,iFrames) = (avgTrace{iCond}(:,:,iFrames) - trace) ./ trace;
    end
end

%% show data for all mixed conditions
for iCond = 1:length(cMagMaps)
    
    [imageOut,areaBounds] = showStack(rot90(avgTrace{iCond},2));
    
    
    outline{iCond} = areaBounds{1};
    set(0,'DefaultFigureWindowStyle','docked')
    figure('name',['stimType = ' num2str(allStimType(iCond))]);
    imagesc(imageOut{1});axis square; colorbar; colormap jet; hold on;
    caxis([0 max(max(imageOut{1}))]);
    plot(smooth(areaBounds{1}(1,:)),smooth(areaBounds{1}(2,:)),'linewidth',2,'color',[1 1 1])
    set(0,'DefaultFigureWindowStyle','normal')
    
end

%% combined figure with whisker and somatosensory stimulus
load('C:\data\WidefieldImager\Animals\mSM21\MixedStim\17-Nov-2016\snapshot_1.mat')
figure(90)
imagesc(rot90(snap,2)); axis square; colormap gray; caxis([0 1000]); freezeColors;
hold on

figure(90); hold on
patch(outline{1}(1,:).*4,outline{1}(2,:).*4,'b','FaceAlpha',0.25); % whisker map
patch(outline{2}(1,:).*4,outline{2}(2,:).*4,'r','FaceAlpha',0.25); % somatosensory map

%% get audiovisual data
load('C:\data\WidefieldImager\Animals\mSM21\MixedStim\17-Nov-2016_1\fftData.mat')
load('C:\data\WidefieldImager\Animals\mSM21\MixedStim\17-Nov-2016_1\avgTrace.mat')

%% normalize
for iCond = 1:length(avgTrace)
    ind = size(avgTrace{iCond},3);
    trace = squeeze(mean(avgTrace{iCond}(:,:,ind/4:ind/2),3));
    for iFrames = 1:size(avgTrace{iCond},3)
        avgTrace{iCond}(:,:,iFrames) = (avgTrace{iCond}(:,:,iFrames) - trace) ./ trace;
    end
end

%% show data for all mixed conditions
set(0,'DefaultFigureWindowStyle','docked')
thresh = [90 70 70];
smth = 2;
clear outline mask

for iCond = 1:length(cMagMaps)
    
    [imageOut,areaBounds] = showStack(rot90(avgTrace{iCond},2));
    
    outline{iCond} = areaBounds;
    set(0,'DefaultFigureWindowStyle','docked')
    figure('name',['stimType = ' num2str(allStimType(iCond))]);
    imagesc(imageOut{1});axis square; colorbar; colormap jet; hold on;
    caxis([0 max(max(imageOut{1}))]);
    plot(smooth(areaBounds{1}(1,:)),smooth(areaBounds{1}(2,:)),'linewidth',2,'color',[1 1 1])
    set(0,'DefaultFigureWindowStyle','normal')
    
end

%% combined figure with audio and visual stimulation
figure(90); hold on
patch(outline{1}{1}(1,:).*4,outline{1}{1}(2,:).*4,'c','FaceAlpha',0.25); % audio map
patch(outline{1}{2}(1,:).*4,outline{1}{2}(2,:).*4,'c','FaceAlpha',0.25); % audio map
patch(outline{2}{1}(1,:).*4,outline{2}{1}(2,:).*4,'y','FaceAlpha',0.25); % vision map

%% get visual data
load('C:\data\WidefieldImager\Animals\mSM21\PhaseMap\17-Nov-2016\plotMap.mat')
load('C:\data\WidefieldImager\Animals\mSM21\PhaseMap\17-Nov-2016\cMagMaps.mat')
temp = imresize(cMagMaps{1},6); %smoothed magnitude map
temp = spatialFilterGaussian(temp,50); %smoothed magnitude map
temp =(temp-min(temp(:)))./(max(temp(:))- min(temp(:))); %normalize between 0 and 1

figure(90); hold on;
vfsIm = imagesc(plotMap); colormap jet; caxis([0.25 0.75]);
set(vfsIm,'AlphaData',temp); axis square

