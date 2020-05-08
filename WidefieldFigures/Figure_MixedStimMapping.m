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
smth = 2;
thresh = [90 90 90 90];
set(0,'DefaultFigureWindowStyle','docked')
for iCond = 1:length(cMagMaps)
    
    temp = smooth2a(cMagMaps{iCond},smth);
    
    figure('name',['stimType = ' num2str(allStimType(iCond))]);
    subplot(2,2,1)
    imagesc(temp);axis square; colorbar; colormap jet; hold on;
    caxis([0 max(max(temp))]);
    freezeColors;
    
    subplot(2,2,2)
    temp = temp;
    temp(:,size(temp,1)/2:end) = NaN;
    [mask{iCond}, outline{iCond}] = getAreas(temp,thresh(iCond),'Area');
    temp(~mask{iCond}) = NaN;
    
    imagesc(smooth2a(temp,smth,smth));axis square; colorbar; colormap jet; hold on;
    caxis([0 max(max(temp))]);
    plot(smooth(outline{iCond}(:,2)),smooth(outline{iCond}(:,1)),'linewidth',2,'color',[1 1 1])
    freezeColors;
    
    subplot(2,2,3)
    temp = smooth2a(cPhaseMaps{iCond},smth);
    imagesc(temp);axis square; colormap hsv; colorbar; hold on;
    caxis([0 max(max(smooth2a(cPhaseMaps{iCond},smth)))]);
    plot(smooth(outline{iCond}(:,2)),smooth(outline{iCond}(:,1)),'linewidth',2,'color',[1 1 1])

    subplot(2,2,4)
    temp(:,size(temp,1)/2:end) = NaN;
    temp(~mask{iCond}) = NaN;
    
    imagesc(temp);axis square; colorbar; hold on;
    caxis([0 max(max(temp))]);
    plot(smooth(outline{iCond}(:,2)),smooth(outline{iCond}(:,1)),'linewidth',2,'color',[1 1 1])
    
    subplot(2,2,1)
    plot(smooth(outline{iCond}(:,2)),smooth(outline{iCond}(:,1)),'linewidth',2,'color',[1 1 1])
        
end
set(0,'DefaultFigureWindowStyle','normal')

%% combined figure with whisker and somatosensory stimulus
load('C:\data\WidefieldImager\Animals\mSM21\MixedStim\17-Nov-2016\snapshot_1.mat')
figure(91)
imagesc(rot90(snap,2)); axis square; colormap gray; caxis([0 1000]); freezeColors;
hold on

figure(91); hold on
patch(outline{1}(:,2).*4,outline{1}(:,1).*4,'b','FaceAlpha',0.25); % whisker map
patch(outline{2}(:,2).*4,outline{2}(:,1).*4,'r','FaceAlpha',0.25); % somatosensory map


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
    
    temp = smooth2a(cMagMaps{iCond},smth);

    figure('name',['stimType = ' num2str(allStimType(iCond))]);
    subplot(2,2,1)
    imagesc(temp); axis square; colorbar; colormap jet; hold on;
    caxis([0 max(max(smooth2a(cMagMaps{iCond},smth)))]); freezeColors;
    
    subplot(2,2,2)
    imagesc(smooth2a(cPhaseMaps{iCond},smth,smth));axis square; colormap hsv; colorbar; hold on;
    caxis([0 max(max(smooth2a(cPhaseMaps{iCond},smth,smth)))]);
    
    subplot(2,2,3)
    temp(:,size(temp,1)/2:end) = NaN;
    temp(smooth2a(cPhaseMaps{iCond},smth) < 0 | smooth2a(cPhaseMaps{iCond},smth) > pi) = NaN;
    [mask{iCond}, outline{iCond}] = getAreas(temp,thresh(iCond),'Area');
    temp(~mask{iCond}) = NaN;
    imagesc(smooth2a(temp,smth,smth));axis square; colorbar; colormap jet; hold on;
    plot(smooth(outline{iCond}(:,2)),smooth(outline{iCond}(:,1)),'linewidth',2,'color',[1 1 1])
    caxis([0 max(max(smooth2a(cMagMaps{iCond},smth)))]);

    subplot(2,2,1)
    plot(smooth(outline{iCond}(:,2)),smooth(outline{iCond}(:,1)),'linewidth',2,'color',[1 1 1])

    subplot(2,2,2)
    plot(smooth(outline{iCond}(:,2)),smooth(outline{iCond}(:,1)),'linewidth',2,'color',[1 1 1])
    
end
set(0,'DefaultFigureWindowStyle','normal')

%% combined figure with whisker and somatosensory stimulus
figure(91); hold on
patch(outline{1}(:,2).*4,outline{1}(:,1).*4,'c','FaceAlpha',0.25); % audio map
patch(outline{2}(:,2).*4,outline{2}(:,1).*4,'y','FaceAlpha',0.25); % vision map
freezeColors;

%% get visual data
load('C:\data\WidefieldImager\Animals\mSM21\PhaseMap\17-Nov-2016\plotMap.mat')
load('C:\data\WidefieldImager\Animals\mSM21\PhaseMap\17-Nov-2016\cMagMaps.mat')
temp = imresize(cMagMaps{1},6); %smoothed magnitude map
temp = spatialFilterGaussian(temp,50); %smoothed magnitude map
temp =(temp-min(temp(:)))./(max(temp(:))- min(temp(:))); %normalize between 0 and 1

figure(91); hold on
vfsIm = imagesc(plotMap); colormap jet; caxis([0.25 0.75]);
set(vfsIm,'AlphaData',temp); axis square

