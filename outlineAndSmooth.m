function outlineSm = outlineAndSmooth(map)

%% Parameters
minArea = 50;
minAreaWidth = 2;

% The window width is what really matters
smoothWindowWidth = 101;
smoothPolyOrder = 2;


%% Produce basic outlines of areas

areas = unique(map(:));
areas(areas == 0) = [];
nAreas = length(areas);

outline = cell(nAreas + 1, 1);
abbrevs = cell(nAreas + 1, 1);
for a = 1:nAreas
  area = imfill(map == areas(a), 'holes');
  outline{a} = bwboundaries(area);
end

% Un-nest cell arrays
outline = vertcat(outline{:});


%% Get rid of junk polygons (tiny areas)
pAreas = NaN(1, length(outline));
pWidth = NaN(1, length(outline));
for p = 1:length(outline)
  pAreas(p) = polyarea(outline{p}(:, 1), outline{p}(:, 2));
  pWidth(p) = min(range(outline{p}(:, 1)), range(outline{p}(:, 2)));
end
outline = outline(pAreas >= minArea & pWidth > minAreaWidth);

%% Smooth
hw = (smoothWindowWidth - 1) / 2;
outlineSm = outline;

% Circularize polygons
for p = 1:length(outlineSm)
  outlineSm{p} = [outlineSm{p}(end-hw+1:end, :); outlineSm{p}; outlineSm{p}(1:hw, :)];
end

% Perform Savitzky-Golay smoothing
for p = 1:length(outlineSm)
  outlineSm{p}(:, 1) = sgolayfilt(outlineSm{p}(:, 1), smoothPolyOrder, smoothWindowWidth);
  outlineSm{p}(:, 2) = sgolayfilt(outlineSm{p}(:, 2), smoothPolyOrder, smoothWindowWidth);
end

% Trim excess points, leave one point of circularization
for p = 1:length(outlineSm)
  outlineSm{p} = outlineSm{p}(hw+1:end-hw+1, :);
end
