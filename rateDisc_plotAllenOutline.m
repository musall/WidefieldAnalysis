function rateDisc_plotAllenOutline(ax,side)
%function to add area outlines from allen CCF to current image

load('allenDorsalMapSM.mat')
if ~exist('ax','var') || isempty(ax)
    ax = gca;
end

cIdx = true(1,length(dorsalMaps.edgeOutlineSplit));
if exist('side','var')
    if strcmpi(side,'L') || strcmpi(side,'R')
        cIdx = ismember(dorsalMaps.sidesSplit,side)';
    end
end

if ishold(ax)
    checker = true;
else
    hold(ax,'on'); checker = false;
end

%check image size
for x = 1 : length(ax.Children)
    if contains(class(ax.Children(x)),'Image')
        lineScale = min((size(dorsalMaps.allenMask) ./ size(ax.Children(x).CData))); %scale lines to match size of the image
    end
end

for x = find(cIdx)
    plot(dorsalMaps.edgeOutlineSplit{x}(:,2)./lineScale, dorsalMaps.edgeOutlineSplit{x}(:,1)./lineScale, 'w', 'LineWidth', 0.1);
end

if ~checker
    hold(ax,'off');
end