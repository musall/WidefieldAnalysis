cPath = 'W:\data\BpodImager\Animals\';
dataOverview = delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);

%% get allen maps
load('allenDorsalMap.mat')
addEdgeOutlinesToDorsalMap(dorsalMaps); %make area figure for orientation
mask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
close(gcf);
f1 = gcf;
outMask = ~poly2mask(dorsalMaps.edgeOutlineSplit{end}(:,2), dorsalMaps.edgeOutlineSplit{end}(:,1),size(mask,1),size(mask,2));
pxPerMM = round(206 / 4);

%% go through recordings and check for alignment
for iAnimals = 1 : length(animals)
    
    if iAnimals == 1 || iAnimals == 13
        h = figure;
    end
    
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep];
    load([fPath 'snapshot_1.mat'])
    load([fPath 'opts2.mat'], 'opts')
    snap = double(snap);
    mask = mask(1:size(snap,1), :); %cut mask to size
    
    % get figure data
    snap = alignAllenTransIm(snap, opts.transParams);   
    snap = snap(:, 1 : size(mask,2));
    cBregma = round(opts.transParams.transBregma);
    aBregma(iAnimals, :) = cBregma;

    % make subplots
    figure(h);
    if iAnimals < 13
        subplot(3, 4, iAnimals)
    else
        subplot(3, 4, iAnimals - 12)
    end
    imagesc(snap);axis image; colormap gray; set(gca,'linewidth',1); axis square
    grid(gca,'on');grid minor;set(gca,'GridColor','w');
    set(gca,'xTick',1:pxPerMM:size(snap,1))
    set(gca,'yTick',1:pxPerMM:size(snap,2))
    hold on;
    
    xVec = cBregma(1)-floor(cBregma(1)/pxPerMM)*pxPerMM:pxPerMM:size(snap,2);
    xLabel = num2str(abs((1:length(xVec))-ceil(cBregma(1)/pxPerMM))');
    set(gca,'xTick',xVec); set(gca,'xTickLabel',xLabel)
    
    yVec = cBregma(2)-floor(cBregma(2)/pxPerMM)*pxPerMM:pxPerMM:size(snap,1);
    yLabel = num2str(abs((1:length(yVec))-ceil(cBregma(2)/pxPerMM))');
    set(gca,'yTick',yVec); set(gca,'yTickLabel',yLabel)
    
    xlabel('ML from Bregma(mm)');
    ylabel('AP from Bregma (mm)');
    
    for x = 1: size(dorsalMaps.edgeOutlineSplit, 1)
        plot(dorsalMaps.edgeOutlineSplit{x}(:,2), dorsalMaps.edgeOutlineSplit{x}(:,1), 'w', 'LineWidth', 0.5);axis image
    end
    plot(round(opts.transParams.transBregma(1)),round(opts.transParams.transBregma(2)), 'o', 'MarkerSize', 10, 'color', 'r', 'linewidth',2)
    drawnow;
    
    % plot bregma into overview figure
    plot(f1.Children, round(opts.transParams.transBregma(1)),round(opts.transParams.transBregma(2)), 'o', 'MarkerSize', 5, 'color', 'r', 'linewidth',2)
end

%% plot average bregma into overview figure
figure(f1);
plot(ceil(mean(aBregma(:,1))), round(mean(aBregma(:,2))), 'o', 'MarkerSize', 5, 'color', 'k', 'linewidth',2);
vline(ceil(mean(aBregma(:,1)))); hline(round(mean(aBregma(:,2))));
disp(mean(aBregma));

%% make another overview figure without area text
figure; hold on
for x = 1: size(dorsalMaps.edgeOutlineSplit, 1)
    plot(dorsalMaps.edgeOutlineSplit{x}(:,2), dorsalMaps.edgeOutlineSplit{x}(:,1), 'k', 'LineWidth', 0.5);axis image
end
plot(aBregma(:,1), aBregma(:,2), 'o', 'MarkerSize', 5, 'color', 'r', 'linewidth',2);
plot(ceil(mean(aBregma(:,1))), round(mean(aBregma(:,2))), 'o', 'MarkerSize', 5, 'color', 'k', 'linewidth',2);
set(gca,'YDir','reverse')
set(gca,'xTick',0:51.5:515); set(gca,'xTickLabel',0 : 1 : 10)
set(gca,'yTick',0:51.5:515); set(gca,'yTickLabel',0 : 1 : 10)
xlabel('ML from Bregma(mm)'); ylabel('AP from Bregma (mm)');

