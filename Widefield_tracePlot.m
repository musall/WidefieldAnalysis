function Widefield_tracePlot(mapData,areaCoord,traceData,mask,trialSegments,cRange,cMap,cIdx,cTitle,xRange,yRange,baseCorrect,grayLine)
% make figure with brain image and individual traces for target coordinates.

if ~exist('cRange', 'var') || isempty(cRange)
    cRange = prctile(abs(mapData(:)),99);
    cRange = [-cRange cRange];
end

if ~exist('cMap', 'var') || isempty(cMap)
    cMap = colormap_blueblackred;
end

if ~exist('trialSegments', 'var')
    trialSegments = [];
end

if ~exist('grayLine', 'var')
    grayLine = false;
end

if ~exist('cIdx', 'var') || isempty(cIdx)
    cIdx = true(1,size(traceData{1},3));
else
    for x = 1:size(traceData,2)
        if length(cIdx) ~= size(traceData{x},3)
          cIdx = true(1,size(traceData{1},3));
          disp('Given index does not match entries of data. Using all entries instead.')
        end
    end
end

if ~exist('cTitle', 'var')
    cTitle = [];
end

if ~exist('xRange', 'var')
    xRange = 1:size(traceData{1},2);
end

if ~exist('yRange', 'var') || isempty(yRange)
    yRange = [-1 2] * 10^-2;
%     yRange = [-1 2];
end

if ~exist('baseCorrect', 'var')
    baseCorrect = false;
end


%%
figure('name','PSTH and traces');
panelIdx = reshape(1:20,5,4)';

subplot(4,5,panelIdx(1:12))
mapImg = imshow(arrayShrink(mapData,mask,'split'),cRange);
traceColor = get(gca,'ColorOrder');
set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
colormap(mapImg.Parent,cMap);hold on; axis image
title(cTitle);

% draw area boundaries, get index and plot traces
traceIdx = cell(1,size(areaCoord,2));
maskIdx = find(~mask);

% reshape data
for iMod = 1:length(traceData)
    if length(cIdx) ~= size(traceData{iMod},3)
        cIdx = true(1,size(traceData{iMod},3));
    end
    for x = 1 : size(areaCoord,2)
        
        areaInfo = regionprops(areaCoord{x} & ~mask, 'Centroid');
        text(areaInfo.Centroid(1)-6,areaInfo.Centroid(2),num2str(x),'color','w','FontSize',15,'parent',mapImg.Parent)
        outline = bwboundaries(areaCoord{x} & ~mask);
        plot(mapImg.Parent,outline{1}(:,2), outline{1}(:,1), 'w', 'LineWidth', 2);
        
        traceIdx{x} = ismember(maskIdx,find(areaCoord{x}));
        regTrace{iMod,x} = squeeze(traceData{iMod}(traceIdx{x},:,cIdx));
        regTrace{iMod,x} = (squeeze(nanmean(regTrace{iMod,x},1)))';
        if baseCorrect  % if traces should be baseline-corrected
%             regTrace{iMod,x} = zscore(regTrace{iMod,x},[],2);
            regTrace{iMod,x} = regTrace{iMod,x} - nanmean(regTrace{iMod,x}(:,1:15),2);
        end
        
        ax = subplot(4,5,panelIdx(x,4:5));hold on
        if isvector(nanmean(regTrace{iMod,x}))
            trace = regTrace{iMod,x};
        else
            trace = nanmean(regTrace{iMod,x});
        end
        if grayLine && iMod == 1
            lines(1,x) = plot(xRange, trace, 'color', [0.5 0.5 0.5], 'linewidth', 2);
        else
            if isvector(trace)
                lines(1,x) = plot(xRange, trace, 'color', traceColor(iMod - grayLine,:), 'linewidth', 2);
            else
                lines(1,x) = stdshade(regTrace{iMod,x},traceColor(iMod - grayLine,:),xRange,0.5);
            end
        end
        xlim([min(xRange) max(xRange)]);
        ylim([min(yRange) max(yRange)]);
        hline(0, '--k'); 
        if ~isempty(trialSegments)
            vline(xRange(trialSegments(2:end)), '--k');
        end
        title(['Area ' num2str(x)]);
        xlabel('time (s)');ylabel('dF/F');
%         axis(ax, 'square');
        
    end
end
end
