function allenMapWithLabels(dorsalMaps, side, useScaled, useEdges)
% allenMapWithLabels(dorsalMaps [, side] [, useScaled] [, useEdges])
% 
% Produce a dorsal-view Allen atlas map, with one hemisphere having the
% areas labeled. Input is dorsalMaps, from the file allenDorsalMap.mat
% (produced by computeAllenDorsalMap).
% 
% By default, the left hemisphere is labeled (side is 'Left' or 'l'). If
% you want the right hemisphere labeled, supply 'Right' or 'r'.
%
% Written by Matt Kaufman, 2018
% 
% To use, please cite:
% Musall S*, Kaufman MT*, Gluf S, Churchland AK (2018). 
% "Movement-related activity dominates cortex during sensory-guided
% decision making." bioRxiv.


%% Optional arguments

if ~exist('side', 'var') || isempty(side) || ~ischar(side)|| strcmpi(side(1), 'l')
  side = 'l';
else
  side = 'r';
end

if ~exist('useScaled', 'var')
  useScaled = 0;
end

if ~exist('useEdges', 'var')
  useEdges = 0;
end


%% Set up map

if ~useScaled
  map = dorsalMaps.dorsalMap;
  dMap = map;
  if useEdges
    dMap = dorsalMaps.edgeMap;
  end
else
  map = dorsalMaps.dorsalMapScaled;
  dMap = map;
  if useEdges
    dMap = dorsalMaps.edgeMapScaled;
  end
end

% Trim
bottom = find(sum(map, 2) > 0, 1, 'last');
map = map(1:bottom, :);
dMap = dMap(1:bottom, :);

% One side only, for centroid finding
if strcmp(side, 'l')
  hMap = map(:, 1:228);
  offset = 0;
else
  hMap = map(:, 229:end);
  offset = 228;
end


%% Identify which areas are visible

uIDs = unique(hMap(:));
uIDs(uIDs == 0) = [];


%% Show map

figure;
imagesc(dMap);

if ~useEdges
  cmap = colorcube(uIDs(end));
  colormap(cmap);
else
  colormap(gray(2));
end

axis equal off;
hold on;


%% Loop through and label each area
for id = uIDs'
  % Find the centroid of the biggest component
  area = (hMap == id);
  rp = regionprops(area, 'centroid', 'area');
  if length(rp) > 1
    [~, biggest] = max([rp.Area]);
    rp = rp(biggest);
  end
  
  % Retrieve the label
  label = getAllenLabel(dorsalMaps.labelTable, id);
  
  % Ensure contrast between label and area color
  if ~useEdges && mean(cmap(id, :)) > 0.5
    color = [0 0 0];
  else
    color = [1 1 1];
  end
  
  text(offset + rp.Centroid(1), rp.Centroid(2), label, 'color', color, 'Horiz', 'center');
end
