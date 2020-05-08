%% get mixedStim data and normalize
cPath = 'C:\data\WidefieldImager\Animals\mSM21\MixedStim\17-Nov-2016\';
load([cPath 'avgTrace.mat'])

for iCond = 1:length(avgTrace)
    ind = size(avgTrace{iCond},3);
    trace = squeeze(mean(avgTrace{iCond}(:,:,ind/4:ind/2),3));
    for iFrames = 1:size(avgTrace{iCond},3)
        avgTrace{iCond}(:,:,iFrames) = (avgTrace{iCond}(:,:,iFrames) - trace) ./ trace;
    end
end

%% combined figure with whisker and somatosensory stimulus
%vessel image
load([cPath 'snapshot_2.mat'])
figure(90)
imshow(imrotate(snap,180)); axis square;
cMap = colormap('gray');
cMap1=cMap; cMap1(:,1:2) = 0; %blue
cMap2=cMap; cMap2(:,[1 3]) = 0; %green
cMap3=cMap; cMap3(:,2:3) = 0; %red
cMap = colormap('gray'); caxis([0 1000]); freezeColors
hold on
% export_fig([cPath 'VesselMap'],'-jpg')

% overlay with mixed conditions
%trunk stim
figure(90)
plotMap = applyFilter2(avgTrace{2}(:,:,18),10,2,2);
plotMap =  imrotate(imresize(plotMap,4),180);
[~,outline{2},~] = getAreas(plotMap,2.5);
trunkMap =(plotMap-min(plotMap(:)))./(max(plotMap(:))- min(plotMap(:))); %normalize between 0 and 1
trunkImg = imagesc(trunkMap); colormap(cMap2); caxis([0 0.75]);
set(trunkImg,'AlphaData',trunkMap); axis square;
freezeColors
% export_fig([cPath 'TrunkProjection'],'-jpg')


% whisker stim
figure(90)
plotMap = applyFilter2(avgTrace{1}(:,:,18),10,2,2);
plotMap =  imrotate(imresize(plotMap,4),180);
[~,outline{1},~] = getAreas(plotMap,4);
barrelMap =(plotMap-min(plotMap(:)))./(max(plotMap(:))- min(plotMap(:))); %normalize between 0 and 1
barrelImg = imagesc(barrelMap); colormap(cMap1); caxis([0 0.75]);
set(barrelImg,'AlphaData',barrelMap); axis square;
freezeColors
% export_fig([cPath 'BarrelProjection'],'-jpg')


%wheel flick
plotMap = applyFilter2(avgTrace{4}(:,:,18),10,2,2);
plotMap =  imrotate(imresize(plotMap,4),180);
[~,outline{3},~] = getAreas(plotMap,2.5);
flickMap =(plotMap-min(plotMap(:)))./(max(plotMap(:))- min(plotMap(:))); %normalize between 0 and 1
figure(90)
% flickImg = imagesc(flickMap); colormap(cMap3); caxis([0 0.75]);
% set(flickImg,'AlphaData',flickMap); axis square;
% freezeColors


%wheel run
plotMap = applyFilter2(avgTrace{3}(:,:,59),10,2,2);
plotMap =  imrotate(imresize(plotMap,4),180);
[~,outline{4},~] = getAreas(plotMap,3);
runMap =(plotMap-min(plotMap(:)))./(max(plotMap(:))- min(plotMap(:))); %normalize between 0 and 1
figure(90)
% runImg = imagesc(runMap); colormap(cMap3); caxis([0 0.75]);
% set(flickImg,'AlphaData',flickMap); axis square;
% freezeColors
% export_fig([cPath 'MixedStimProjection'],'-jpg')
% 
% plot(smooth(outline{1}(:,2),10),smooth(outline{1}(:,1),10),'b','linewidth',2); %whisker
% plot(smooth(outline{2}(:,2),10),smooth(outline{2}(:,1),10),'g','linewidth',2); %trunk
% plot(smooth(outline{3}(:,2),10),smooth(outline{3}(:,1),10),'k','linewidth',2); %wheel run
% plot(smooth(outline{4}(:,2),10),smooth(outline{4}(:,1),10),'r','linewidth',2); %wheel flick

% export_fig([cPath 'MixedStimProjection'],'-pdf')

%% get visual data
load('C:\data\WidefieldImager\Animals\mSM21\PhaseMap\17-Nov-2016\plotPhaseMap.mat')
load('C:\data\WidefieldImager\Animals\mSM21\PhaseMap\17-Nov-2016\cMagMaps.mat')
temp = imresize(cMagMaps{1},4); %smoothed magnitude map
temp = spatialFilterGaussian(temp,50); %smoothed magnitude map
temp =(temp-min(temp(:)))./(max(temp(:))- min(temp(:))); %normalize between 0 and 1

figure(90); hold on;
vfsIm = imagesc(plotPhaseMap); colormap jet; caxis([0.25 0.75]);
set(vfsIm,'AlphaData',temp); axis square
% export_fig([cPath 'MixedStimProjection_VFS'],'-jpg')

%% show individual responses with outline
figure(50);
subplot(2,2,1)
imshow(barrelMap); hold on; title('Whisker'); caxis([0 0.85])
plot(smooth(outline{1}(:,2),10),smooth(outline{1}(:,1),10),'k:','linewidth',1);

subplot(2,2,2)
imshow(trunkMap); hold on; title('Trunk'); caxis([0 0.85])
plot(smooth(outline{2}(:,2),10),smooth(outline{2}(:,1),10),'k:','linewidth',1);

subplot(2,2,3)
imshow(flickMap); hold on; title('WheelFlick'); caxis([0 0.85])
plot(smooth(outline{3}(:,2),10),smooth(outline{3}(:,1),10),'k:','linewidth',1);

subplot(2,2,4)
imshow(runMap); hold on; title('WheelRun'); caxis([0 0.85])
plot(smooth(outline{4}(:,2),10),smooth(outline{4}(:,1),10),'k:','linewidth',1);

% export_fig([cPath 'MixedStimProjection_subplots'],'-pdf')

%% 
close(51)
close(90);

%% get mixedStim data and normalize
cPath = 'C:\data\WidefieldImager\Animals\mSM21\MixedStim\17-Nov-2016_1\';
load([cPath 'avgTrace.mat'])

for iCond = 1:length(avgTrace)
    ind = size(avgTrace{iCond},3);
    trace = squeeze(mean(avgTrace{iCond}(:,:,ind/4:ind/2),3));
    for iFrames = 1:size(avgTrace{iCond},3)
        avgTrace{iCond}(:,:,iFrames) = (avgTrace{iCond}(:,:,iFrames) - trace) ./ trace;
    end
end

%% combined figure with whisker and somatosensory stimulus
%vessel image
load('C:\data\WidefieldImager\Animals\mSM21\MixedStim\17-Nov-2016\snapshot_1.mat')
figure(90)
imshow(imrotate(snap,180)); axis square;
cMap = colormap('gray');
cMap1=cMap; cMap1(:,1:2) = 0; %blue
cMap2=cMap; cMap2(:,[1 3]) = 0; %green
cMap3=cMap; cMap3(:,2:3) = 0; %red
cMap = colormap('gray'); caxis([0 1000]); freezeColors
hold on

%% overlay with mixed conditions
%audio stim
figure(90)
plotMap = applyFilter2(avgTrace{1}(:,:,18),10,2,2);
plotMap =  imrotate(imresize(plotMap,4),180);
temp = plotMap; temp(:,size(temp,1)/2:end) = 0; %exclude right HS
[~,outline{1},~] = getAreas(temp,2.5); clear temp
audioMap =(plotMap-min(plotMap(:)))./(max(plotMap(:))- min(plotMap(:))); %normalize between 0 and 1
audioImg = imagesc(audioMap); colormap(cMap3); caxis([0 1]);
set(audioImg,'AlphaData',audioMap); axis square;
freezeColors
% export_fig([cPath 'AudioProjection'],'-jpg')

plotMap = applyFilter2(avgTrace{2}(:,:,17),10,2,2);
plotMap =  imrotate(imresize(plotMap,4),180);
temp = plotMap; temp(:,size(temp,1)/2:end) = 0; %exclude right HS
[~,outline{2},~] = getAreas(temp,1); clear temp
visionMap =(plotMap-min(plotMap(:)))./(max(plotMap(:))- min(plotMap(:))); %normalize between 0 and 1
visionImg = imagesc(visionMap); colormap(cMap1); caxis([0 .75]);
set(visionImg,'AlphaData',visionMap); axis square;
freezeColors
% export_fig([cPath 'VisionProjection'],'-jpg')


% plot(smooth(outline{1}(:,2),10),smooth(outline{1}(:,1),10),'r','linewidth',2); %audio
% plot(smooth(outline{2}(:,2),10),smooth(outline{2}(:,1),10),'b','linewidth',2); %vision

%% get visual data
load('C:\data\WidefieldImager\Animals\mSM21\PhaseMap\17-Nov-2016\plotPhaseMap.mat')
load('C:\data\WidefieldImager\Animals\mSM21\PhaseMap\17-Nov-2016\cMagMaps.mat')
temp = imresize(cMagMaps{1},4); %smoothed magnitude map
temp = spatialFilterGaussian(temp,50); %smoothed magnitude map
temp =(temp-min(temp(:)))./(max(temp(:))- min(temp(:))); %normalize between 0 and 1

figure(90); hold on;
vfsIm = imagesc(plotPhaseMap); colormap jet; caxis([0.25 0.75]);
set(vfsIm,'AlphaData',temp); axis square
% export_fig([cPath 'AudioVisualProjection_VFS'],'-jpg')
% export_fig([cPath 'AllCombinedProjection_VFS'],'-jpg')

%% show individual responses with outline
figure(51);
subplot(2,2,1)
imshow(audioMap); hold on; title('Whisker'); caxis([0 0.85])
plot(smooth(outline{1}(:,2),10),smooth(outline{1}(:,1),10),'k:','linewidth',1);

subplot(2,2,2)
imshow(visionMap); hold on; title('Vision'); caxis([0 0.5])
plot(smooth(outline{2}(:,2),10),smooth(outline{2}(:,1),10),'w:','linewidth',1);

% export_fig([cPath 'MixedAudioVisualProjection_subplots'],'-pdf')
