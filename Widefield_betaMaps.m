% function BpodImager_betaMaps(cPath)
cPath = 'U:\smusall\BpodImager\Animals\';
% cPath = '\\CHURCHLANDNAS\homes\DOMAIN=CSHL\smusall\NNpaper_Dataset\Widefield\'; %path to local NAS
[dataOverview, motorLabels, sensorLabels, cogLabels] = delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);
taskLabels = [sensorLabels cogLabels]; %combine sensory and cognitive labels, so sensory is really 'non-motor' or 'task' here
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
audExp = ismember(dataOverview(:,2),'Audio'); %index for auditory experts
videoLength = 90; %nr of frames per beta kernel
% animals = animals(audExp); %only visual animals
% recs = recs(audExp); %only visual animals

%%
for iAnimals = 1 : length(animals)
    tic
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep]; %path to current dataset
%     disp(fPath); %current data path

    % get some data
    load([fPath 'dimBeta.mat'])
    load([fPath 'regData.mat'], 'recIdx', 'idx', 'recLabels')
    load([fPath 'Vc.mat'], 'U')
    load([fPath 'opts2.mat'], 'opts')
            
    if iAnimals == 1
        load('allenDorsalMap.mat')
        load([fPath 'snapshot_1.mat'])
        mask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
        [x1, y1] = size(mask);
        addEdgeOutlinesToDorsalMap(dorsalMaps); %make area figure for orientation
        mask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
        snap = alignAllenTransIm(single(snap),opts.transParams);
        [x2, y2] = size(snap);
        mask = mask(1:min([x1 x2]), 1:min([y1 y2])); %cut mask to size
        mask = mask(1:min([x1 x2]), 1:min([y1 y2])); %cut mask to size
        rightHs = find(ismember(dorsalMaps.sidesSplit,'R')); %index for labels on the left HS
        allVidVar = NaN(length(animals),sum(recIdx(~idx) == find(strcmpi(recLabels, 'bhvVideo')))-1, 'single'); %pre-allocate video eigenvalue matrix
    end
    
    U = U(:, :, 1:size(dimBeta,2));
    U = alignAllenTransIm(U, opts.transParams);
    U = arrayShrink(U(1:size(mask,1), 1:size(mask,2),:), mask);
    
    if iAnimals ==1 
        allU = NaN(size(U,1), size(U,2), length(animals), 'single');
    end
    allU(:,:,iAnimals) = U; clear U
    
    % get dimensions for different betas
    for iRegs = unique(recIdx)
        allRegs{iAnimals, iRegs} = recLabels{iRegs};
        betaData{iAnimals, iRegs} = NaN(sum(recIdx == iRegs), size(dimBeta,2), 'single'); %pre-allocate for all weights
        betaData{iAnimals, iRegs}(~idx(recIdx == iRegs),:) = dimBeta(recIdx(~idx) == iRegs, :);  %pick up used regressor weights
    end
    
    if rem(iAnimals, round(length(animals)/5)) == 0
        fprintf(1, 'Done. Current recording is %d out of %d\n', iAnimals, length(animals));
        toc;
    end
    
    % run pca on video weights
    vidBeta = dimBeta(recIdx(~idx) == find(strcmpi(recLabels, 'bhvVideo')), :); %video weights
    [~,~,c] = pca(vidBeta);
    allVidVar(iAnimals,:) = 1-c./sum(c);
end

%% Compute beta maps - averaged over all weights
% vidLabels = {'lVisStim' 'rVisStim' 'lAudStim' 'rAudStim' 'lGrab' 'rGrab' 'nose' 'lLick' 'rLick' 'whisk' 'piezo'};
vidLabels = {'rVisStim' 'rGrab' 'nose' 'whisk' 'piezo' 'fastPupil'};
betaVideos = NaN(size(mask,1), size(mask,2), videoLength, length(vidLabels), 'single');

for iRegs = 1 : length(recLabels)
    for iAnimals = 1 : length(animals)
        
        if ismember(recLabels{iRegs}, vidLabels)
            cLength = size(betaData{iAnimals, iRegs},1);
            if ismember(recLabels{iRegs},{'nose'})
                cRange = ceil((cLength-1)/2) + 2 : cLength; %high amp events
%             elseif ismember(recLabels{iRegs},{'whisk' 'piezo' })
            elseif ismember(recLabels{iRegs},{'whisk'})
                cRange = 2 : ceil((cLength-1)/2) + 1; %low amp events
            else
                cRange = 1 : cLength;
            end
            if length(cRange) > videoLength
                cRange = cRange(1 : videoLength);
            end
            
            a = betaVideos(:,:,1:length(cRange),ismember(vidLabels, recLabels{iRegs})) .* (iAnimals-1);
            b = arrayShrink(allU(:,:,iAnimals) * betaData{iAnimals, iRegs}(cRange,:)', mask, 'split');

            if any(abs(b(:)) > nanstd(b(:))*8) %check for weird edges and erode a little if needed
                temp = ~imerode(~isnan(b(:,:,1)),strel('square', 3));
                b = arrayShrink(b, temp);
                b = arrayShrink(b, temp,'split');
            end
            c = a + b;
            c(isnan(c)) = a(isnan(c));
            c(isnan(c)) = b(isnan(c));
            betaVideos(:,:,1:length(cRange),ismember(vidLabels, recLabels{iRegs})) = c ./ iAnimals; 
            clear a b c
        end
        
    end
end

%% show scree plot for video weights
set(0, 'DefaultFigureRenderer', 'painters');
figure
cData = 1-cumsum(1-allVidVar(:,1:10),2); 
plot(cData','color',[.5 .5 .5]);hold on; axis square
errorbar(mean(cData), sem(cData), '-k', 'linewidth', 2); axis square
plot(mean(cData), 'ok', 'linewidth', 2, 'MarkerSize', 6, 'MarkerFaceColor','w');
hline(0.1); ylabel('Explained variance'); xlabel('Components');
title('Video weight cummulative variance');
xlim([0.5 10.5]);

%% show handles and visual stim beta
sRate = 30;
cParams = {'rVisStim' 'rGrab' 'nose' 'whisk'}; %variables to show
mapTimes = {[0.2 1],[0 0.333]+0.5,[0 0.333]+0.5,[0.1 0.333]+0.5}; %timepoints to show. should be 3 per variable.
areaIdx = {'VISp' 'SSp-ul' 'MOB' 'SSp-bfd'}; % V1 - HL - M2
% get allen area masks
for x = 1:length(areaIdx)
    ind = ismember(dorsalMaps.labelsSplit,areaIdx{x}) & ismember(dorsalMaps.sidesSplit,'L');
    areaCoord{x} = poly2mask(dorsalMaps.edgeOutlineSplit{ind}(:,2), dorsalMaps.edgeOutlineSplit{ind}(:,1),size(mask,1),size(mask,2));
end

figure
cRange = 0.01;
for x = 1:2
subplot(3,3,x);
cMap = betaVideos(:,:,round(mapTimes{1}(x)*sRate),ismember(vidLabels,cParams{1}));
mapImg = imshow(cMap,[-cRange cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make 0s transparent.
colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
title([cParams{1} '; ' num2str(mapTimes{1}(x)) 's']); colorbar
hold(mapImg.Parent, 'on');
for xx = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{xx}(:,2), dorsalMaps.edgeOutlineSplit{xx}(:,1), 'w', 'LineWidth', 0.2);axis image
end
end
cTrace = nanmean(arrayShrink(betaVideos(:,:,:,ismember(vidLabels,cParams{1})), ~areaCoord{1}, 'merge'),1);
subplot(3,3,3); plot(0 : 1/sRate : 81/sRate,cTrace(1:82),'k','linewidth',2);
xlabel('Time(s)'); ylabel('weight');xlim([-0.1 2.75]);
ylim([-0.005 0.02]);

y = 2;
cRange = 0.0015;
for x = 1:2
subplot(3,3,(y-1)*3+x);
cMap = betaVideos(:,:,round(mapTimes{y}(x)*sRate),ismember(vidLabels,cParams{y}));
mapImg = imshow(cMap,[-cRange cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make 0s transparent.
colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
title([cParams{y} '; ' num2str(mapTimes{y}(x)-0.5) 's']); colorbar
for xx = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{xx}(:,2), dorsalMaps.edgeOutlineSplit{xx}(:,1), 'w', 'LineWidth', 0.2);axis image
end
end
cTrace = nanmean(arrayShrink(betaVideos(:,:,:,ismember(vidLabels,cParams{2})), ~areaCoord{2}, 'merge'),1);
subplot(3,3,6); plot(-0.5 : 1/sRate : 66/sRate,cTrace(1:82),'k','linewidth',2);
xlabel('Time(s)'); ylabel('weight');xlim([-0.55 2]);
ylim([-0.001 0.0025]); 
vline(0);

y = 3;
cRange = 0.0005;
for x = 1:2
subplot(3,3,(y-1)*3+x);
cMap = betaVideos(:,:,round(mapTimes{y}(x)*sRate),ismember(vidLabels,cParams{y}));
mapImg = imshow(cMap,[-cRange cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make 0s transparent.
colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
title([cParams{y} '; ' num2str(mapTimes{y}(x)-0.5) 's']); colorbar
for xx = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{xx}(:,2), dorsalMaps.edgeOutlineSplit{xx}(:,1), 'w', 'LineWidth', 0.2);axis image
end
end
cTrace = nanmean(arrayShrink(betaVideos(:,:,:,ismember(vidLabels,cParams{3})), ~areaCoord{3}, 'merge'),1);
subplot(3,3,9); plot(-0.5 : 1/sRate : 66/sRate,cTrace(1:82),'k','linewidth',2);
xlabel('Time(s)'); ylabel('weight');xlim([-0.55 2]);
ylim([-0.0005 0.00075]);
vline(0);


%% pupil beta map
figure;
cMap = betaVideos(:,:,1,ismember(vidLabels,'fastPupil'));
cRange = 0.0004;
mapImg = imshow(cMap,[-cRange cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make 0s transparent.
colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
title('Pupil'); colorbar


%% whisker filter example figure (use regData)
y = 4;
cParams = {'rVisStim' 'rGrab' 'nose' 'whisk'}; %variables to show
whiskTimes = [0.14 0.6 0.86]; %timepoints to show. should be 3 per variable.
whiskAreas = {'SSp-bfd' 'RSPd'}; %BCX - RS
% get allen area masks
for x = 1:length(whiskAreas)
    ind = ismember(dorsalMaps.labelsSplit,whiskAreas{x}) & ismember(dorsalMaps.sidesSplit,'L');
    whiskCoord{x} = poly2mask(dorsalMaps.edgeOutlineSplit{ind}(:,2), dorsalMaps.edgeOutlineSplit{ind}(:,1),size(mask,1),size(mask,2));
end

figure
cRange = 0.0003;
for x = 1:3
subplot(2,2,x);
cMap = betaVideos(:,:,round(whiskTimes(x)*sRate),ismember(vidLabels,cParams{y}));
mapImg = imshow(cMap,[-cRange cRange]);
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make 0s transparent.
colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
title([cParams{y} '; ' num2str(whiskTimes(x)-0.5) 's']); colorbar
for xx = 1: length(dorsalMaps.edgeOutlineSplit)
    plot(dorsalMaps.edgeOutlineSplit{xx}(:,2), dorsalMaps.edgeOutlineSplit{xx}(:,1), 'w', 'LineWidth', 0.2);axis image
end
end
cTrace = nanmean(arrayShrink(betaVideos(:,:,:,ismember(vidLabels,cParams{4})), ~whiskCoord{1}, 'merge'),1);
subplot(2,2,4); plot(-0.5 : 1/sRate : 66/sRate,cTrace(1:82),'k','linewidth',2); hold on;
cTrace = nanmean(arrayShrink(betaVideos(:,:,:,ismember(vidLabels,cParams{4})), ~whiskCoord{2}, 'merge'),1);
subplot(2,2,4); plot(-0.5 : 1/sRate : 66/sRate,cTrace(1:82),'r','linewidth',2); hold on;
xlabel('Time(s)'); ylabel('weight');xlim([-0.55 2]);
ylim([-0.025 0.0275]); vline(0); 
vline(whiskTimes-0.5)


%% make beta movies
moviePath = 'X:\smusall\betaMovies';
for iRegs = 1 : length(vidLabels)
    cMovie = double(betaVideos(:,:,1:75,iRegs));
    disp(vidLabels(iRegs))
    compareMovie(cMovie,'label',vidLabels(iRegs))
end

%% analog regressor example
p1 = 25; % time of first pulse
p2 = 700; %second pulse time
shift = 200; %time shift for shifted example
x = 1:1000; %filter length

gauss1 = normpdf(x,300,50); 
gauss1 = gauss1 ./ max(gauss1);

gam1 = gampdf(x,2,100); 
gam1 = gam1 ./ max(gam1);

bp1 = normpdf(x,300,50) - normpdf(x,450,100); 
bp1 = bp1 ./ max(bp1);

trace1 = zeros(1,1000);
trace2 = zeros(1,1000);
trace1(p1) = 1; 
trace2(p2) = 1;

sTrace1 = zeros(1,1000);
sTrace2 = zeros(1,1000);
sTrace1(p1+shift) = 1; 
sTrace2(p2+shift) = 1;

a = conv(trace1,gauss1) + conv(trace2,gauss1)./3;
b = conv(sTrace1,gauss1) + conv(sTrace2,gauss1)./3;
c = conv(sTrace1,gam1) + conv(sTrace2,gam1)./3;
d = conv(trace1,bp1) + conv(trace2,bp1)./3;

subplot(1,5,1);
plot(a); axis square
subplot(1,5,2);
plot(a); hold on; plot(a*2); hold off; ylim([0 2.5]); axis square
subplot(1,5,3);
plot(a); hold on; plot(b); hold off; ylim([0 1.25]); axis square
subplot(1,5,4);
plot(a); hold on; plot(c); hold off; ylim([0 1.25]); axis square
subplot(1,5,5);
plot(a); hold on; plot(d); hold off; ylim([-1.5 1.5]); axis square


%% summary figure animation
rowLabels{1} = {'water' 'prevMod' 'prevChoice' 'prevReward'};
rowLabels{2} = {'choice' 'visReward' 'audReward' 'time'};
rowLabels{3} = {'lVisStim' 'rVisStim' 'lAudStim' 'rAudStim'};


%% Make task figure

rowLabels{1} = {'water' 'prevMod' 'prevChoice' 'prevReward'};
rowLabels{2} = {'choice' 'visReward' 'audReward' 'time'};
rowLabels{3} = {'lVisStim' 'rVisStim' 'lAudStim' 'rAudStim'};

figure
Cnt = 0;
for iRows = 1 : length(rowLabels)
    for iMaps = 1 : length(rowLabels{iRows})
        cMap = betaMaps(:, :, strcmp(recLabels, rowLabels{iRows}{iMaps}));
        allScale(iRows, iMaps) =  max(abs(prctile(cMap(:),[1 99])));
    end
    
    for iMaps = 1 : length(rowLabels{iRows})
        Cnt = Cnt + 1;
        subplot(length(rowLabels), 4 ,Cnt)
        
        cMap = betaMaps(:, :, strcmp(recLabels, rowLabels{iRows}{iMaps}));
        allScale(iRows, iMaps) =  max(abs(prctile(cMap(:),[1 99])));
        cRange = mean(allScale(iRows,:));
        mapImg = imshow(cMap,[-cRange cRange]);
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
        title(rowLabels{iRows}{iMaps});
        drawnow;
        
        hold(mapImg.Parent, 'on');
        for x = 1: length(rightHs)
            plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
        end
        
    end
end


%% Make movement figure

rowLabels{1} = {'whisk' 'nose' 'pupil' 'lLick'};
rowLabels{2} = {'body' 'piezo' 'lGrab' 'rGrab'};

figure
Cnt = 0;
for iRows = 1 : length(rowLabels)
    for iMaps = 1 : length(rowLabels{iRows})
        cMap = betaMaps(:, :, strcmp(recLabels, rowLabels{iRows}{iMaps}));
        allScale(iRows, iMaps) =  max(abs(prctile(cMap(:),[1 99])));
    end
    
    for iMaps = 1 : length(rowLabels{iRows})
        Cnt = Cnt + 1;
        subplot(length(rowLabels), 4 ,Cnt)
        
        cMap = betaMaps(:, :, strcmp(recLabels, rowLabels{iRows}{iMaps}));
        allScale(iRows, iMaps) =  max(abs(prctile(cMap(:),[1 99])));
        cRange = mean(allScale(iRows,:));
        mapImg = imshow(cMap,[-cRange cRange]);
        set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
        colormap(mapImg.Parent,colormap_blueblackred(256));hold on;
        title(rowLabels{iRows}{iMaps});
        drawnow;
        
        hold(mapImg.Parent, 'on');
        for x = 1: length(rightHs)
            plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
        end
        
    end
end