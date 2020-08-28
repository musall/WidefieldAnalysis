function alignBrainToAllen(varargin)
% alignBrainToAllen(WFIm)
% alignBrainToAllen(animal, theDate)
% 
% GUI to align brains to the Allen atlas.
% 
% Inputs:
% Option 1: WFIm should be some useful image from the WF, such as blueAvg.
%           Once alignment is done, you can send the resulting parameters
%           to the Matlab workspace.
% Option 2: animal and theDate specify which dataset to load. Automatically
%           selects blueAvg.mat as the image. With this option, at the end
%           you'll get a button to save the results to opts2.mat (or you
%           can send the result to the workspace, as with option 1).


%% Parse inputs

haveOpts2 = 0;

if nargin == 1
  dat.WFIm = varargin{1};
  dat.havePath = 0;
  
elseif nargin == 2
  dat.havePath = 1;
  
  if ismac
    basePath = '/Volumes/churchland_hpc_home/space_managed_data/BpodImager/Animals';
  elseif ispc
    basePath = 'U:\smusall\BpodImager\Animals';
%     basePath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals';
  end
  
  if ~exist(basePath, 'dir')
    error('Could not access path; server probably not mounted. Tried: %s', basePath);
  end
  
  dat.path = fullfile(basePath, varargin{1}, 'SpatialDisc', varargin{2});
  
  if exist(fullfile(dat.path, 'blueAvg.mat'), 'file')
    load(fullfile(dat.path, 'blueAvg.mat'));
  else
    error('Could not find blueAvg.mat');
  end
  
  dat.WFIm = blueAvg;
  
  if exist(fullfile(dat.path, 'opts.mat'), 'file')
    % If we get here, we don't need to throw an error but may still not
    % have opts2.mat
    if exist(fullfile(dat.path, 'opts2.mat'), 'file')
      load(fullfile(dat.path, 'opts2.mat'));
      haveOpts2 = 1;
    end
  else
    error('Could not find opts.mat');
  end

else
  error('Supply either an image or an animal name and dataset date');
end


%% Parameters

% Points to click, in Allen coordinates relative to our trimmed brain image
dat.refPts = [
  148 79;    % Base of L OB
  228.5 76   % Center base of OBs
  308 79     % Base of R OB
  228.5 335  % Base of RS
  ];

dat.bregmaRef = [228.5 190];

dat.allenCrossSize = 15;
dat.allenCrossLineWidth = 3;
dat.crossColors = {'g', 'r', 'c', 'w'};
dat.bregmaColor = 'y';

dat.selectedCrossSize = 12;


dat.edgeColor = [1 1 1];



%% Load Allen data maps

loadvar = load('allenDorsalMap', 'dorsalMaps');
dat.aMaps = loadvar.dorsalMaps;
dat.nAreas = max(dat.aMaps.dorsalMap(:));


%% Set up figure

hf = figure;
hf.Position = [100 250 1200 500];
hf.NumberTitle = 'off';

hf.Name = 'Brain Aligner';
if dat.havePath
  hf.Name = sprintf('Brain Aligner: %s %s', varargin{1}, varargin{2});
end

drawnow;
hf.Pointer = 'crosshair';


%% Set up reference axes, plot Allen brain

% Axes
dat.axAllen = axes(hf, 'Position', [0.05 0.05 0.35 0.9]);

% Plot brain
imagesc(dat.aMaps.dorsalMap);
colormap(colorcube(dat.nAreas));
axis equal off;
hold on;


%% Set up widefield image axes

dat.axWF = axes(hf, 'Position', [0.5 0.05 0.35 0.9]);

maxBright = prctile(dat.WFIm(:), 99);
dat.WFIm(dat.WFIm > maxBright) = maxBright;

hIm = imagesc(dat.WFIm);
colormap(dat.axWF, gray(256));
axis equal off;
hold on;

% Make brain clickable
hIm.ButtonDownFcn = @advanceSelection;


%% "Allow scaling" checkbox

dat.chkScale = uicontrol('Style', 'checkbox', ...
  'String', 'Allow scaling', ...
  'Units', 'Normalized', 'Position', [0.87 0.5 0.1 0.04], 'Value', 1);


%% "Done" button

dat.btnDone = uicontrol('Style', 'pushbutton', ...
  'String', 'Done', ...
  'Units', 'Normalized', 'Position', [0.87 0.05 0.1 0.04], ...
  'Callback', @performAlignment, 'Enable', 'off');


%% Set up data structures to hold selections and reference crosses

if haveOpts2
  % Use previous selections
  dat.selectedPts = opts.transParams.selectedPts';
  dat.bregma = opts.transParams.selectedBregma;
  
  % Draw objects
  for ri = 1:size(dat.refPts, 1)
    % Reference crosses
    plot(dat.axAllen, dat.refPts(ri, 1), dat.refPts(ri, 2), '+', ...
      'MarkerSize', dat.allenCrossSize, 'LineWidth', dat.allenCrossLineWidth, ...
      'color', dat.crossColors{ri});
    
    % Selected points
    dat.WFCrossHs(ri) = plot(dat.axWF, dat.selectedPts(ri, 1), dat.selectedPts(ri, 2), '+', ...
      'MarkerSize', dat.selectedCrossSize, 'color', dat.crossColors{ri});
    
    % Make cross draggable
    draggableActPlusUp(dat.WFCrossHs(ri), [], @moveCross);
  end
  
  % Reference Bregma
  plot(dat.axAllen, dat.bregmaRef(1), dat.bregmaRef(2), '+', ...
    'MarkerSize', dat.allenCrossSize, 'LineWidth', dat.allenCrossLineWidth, ...
    'color', dat.bregmaColor);
  plot(dat.axAllen, dat.bregmaRef(1), dat.bregmaRef(2), 'o', ...
    'MarkerSize', dat.allenCrossSize - 4, 'LineWidth', dat.allenCrossLineWidth, ...
    'color', dat.bregmaColor);
  
  % Selected Bregma
  hc = plot(dat.axWF, dat.bregma(1), dat.bregma(2), 'x', ...
    'MarkerSize', dat.selectedCrossSize, 'color', dat.bregmaColor);
  
  % Make cross draggable
  draggableActPlusUp(hc, [], @moveBregma);
  
  
  dat.refI = size(dat.refPts, 1) + 1;
  dat.btnDone.Enable = 'on';
  
  if opts.transParams.scaleConst == 1
    dat.chkScale.Value = 0;
  else
    dat.chkScale.Value = 1;
  end
  
  hf.Pointer = 'arrow';
  
  hf.UserData = dat;

  title(dat.axWF, 'Drag crosses to revise');

else
  % Set up first selection
  
  dat.selectedPts = NaN(size(dat.refPts));
  dat.WFCrossHs = [];
  dat.refI = 1;
  hf.UserData = dat;
  
  % Draw first target cross
  drawCrossOnAllen(hf);
  
  title(dat.axWF, 'Click to place a ref point');
end




function drawCrossOnAllen(hf)

%% Retrieve all the data from the figure

dat = hf.UserData;

if dat.refI <= size(dat.refPts, 1)
  % Normal reference point
  plot(dat.axAllen, dat.refPts(dat.refI, 1), dat.refPts(dat.refI, 2), '+', ...
    'MarkerSize', dat.allenCrossSize, 'LineWidth', dat.allenCrossLineWidth, ...
    'color', dat.crossColors{dat.refI});
else
  % Bregma
  plot(dat.axAllen, dat.bregmaRef(1), dat.bregmaRef(2), '+', ...
    'MarkerSize', dat.allenCrossSize, 'LineWidth', dat.allenCrossLineWidth, ...
    'color', dat.bregmaColor);  
  plot(dat.axAllen, dat.bregmaRef(1), dat.bregmaRef(2), 'o', ...
    'MarkerSize', dat.allenCrossSize - 4, 'LineWidth', dat.allenCrossLineWidth, ...
    'color', dat.bregmaColor);  
end

hf.UserData = dat;


function advanceSelection(src, ev)

%% Retrieve all the data from the figure

hf = gcbf;
dat = hf.UserData;


if strcmp(hf.SelectionType, 'normal')
  %% Left click, advance
  
  % If we've made as many selections as we can (including Bregma), ignore click
  if dat.refI > size(dat.refPts, 1) + 1
    return;
  end
  
  
  %% Find click
  
  % Current mouse coordinates
  pt = dat.axWF.CurrentPoint(1, 1:2);
  
  if dat.refI <= size(dat.refPts, 1)
    % Normal reference point
    dat.selectedPts(dat.refI, :) = pt;
    dat.WFCrossHs(dat.refI) = plot(dat.axWF, pt(1), pt(2), '+', ...
      'MarkerSize', dat.selectedCrossSize, 'color', dat.crossColors{dat.refI});
    
    % Make cross draggable
    draggableActPlusUp(dat.WFCrossHs(dat.refI), [], @moveCross);
  else
    % Bregma
    dat.bregma = pt;
    hc = plot(dat.axWF, pt(1), pt(2), 'x', ...
      'MarkerSize', dat.selectedCrossSize, 'color', dat.bregmaColor);
    
    % Make cross draggable
    draggableActPlusUp(hc, [], @moveBregma);
  end
  
  
  %% Advance
  
  dat.refI = dat.refI + 1;
  
  hf.UserData = dat;
  
  if dat.refI > size(dat.refPts, 1) + 1
    hf.Pointer = 'arrow';
    dat.btnDone.Enable = 'on';
    title(dat.axWF, 'Drag markers to revise');
  elseif dat.refI == size(dat.refPts, 1) + 1
    drawCrossOnAllen(hf);
    title(dat.axWF, 'Click to estimate Bregma (ignore the ''V'' in the suture)');
  else
    drawCrossOnAllen(hf);
    title(dat.axWF, 'Click to place a ref point or drag a cross to revise');
  end
end


function performAlignment(src, ev)

hf = gcbf;
dat = hf.UserData;

%% Rescale reference points to our resolution

scaledRefs = (dat.refPts - 1) * dat.aMaps.allenPixelSize / dat.aMaps.desiredPixelSize + 1;


%% Perform the alignment itself, in 2 coordinate frames

% Image coordinates:
% Use the Center of the image, because that's how imrotate works. The +0.5
% is because Matlab indexing starts at 1.
c = fliplr(size(dat.WFIm)) / 2 + 0.5;
scaledRefsC = scaledRefs - c;
selectedPtsC = dat.selectedPts - c;

% Kabsch algorithm to align the points with a rotation and translation
[RC, tC] = rigidAlignPts(selectedPtsC', scaledRefsC', dat.chkScale.Value);   % For image
[R, t] = rigidAlignPts(dat.selectedPts', scaledRefs', dat.chkScale.Value);   % For alignment points

% Retrieve the angle from the rotation matrix, rotate the image
transParams.scaleConst = sqrt(R(1, 1) ^ 2 + R(1, 2) ^ 2);
ang = acosd(RC(1) / transParams.scaleConst);
ang = ang * sign(RC(1, 2));

% Pack the outputs
transParams.angleD = ang;
transParams.R = R;
transParams.t = t;
transParams.RC = RC;
transParams.tC = tC;

WF = alignAllenTransIm(dat.WFIm, transParams);


%% Set up figure

hf = figure;
hf.Position = [100 150 1200 500];


%% First panel: Aligned image with selected points and reference points

% Annoyingly, we need to know the colormap resolution used inside
% showAllenEdgesOnWF, because it will override the colormap
cRes = 255;

subplot(1, 3, 1);
image(cRes * (WF - min(WF(:))) / range(WF(:)));
% colormap(gray);
axis equal off;

% Crosses for reference points, X's for selected points. This lets us see
% how well they aligned
hold on;
plot(scaledRefs(2, 1) * [1 1], [1 size(WF, 1)], 'w-');
selPts = (R * dat.selectedPts' + t)';
for p = 1:size(dat.refPts, 1)
  plot(scaledRefs(p, 1), scaledRefs(p, 2), '+', ...
    'MarkerSize', 12, 'LineWidth', 1, 'color', dat.crossColors{p});
  plot(selPts(p, 1), selPts(p, 2), 'x', ...
    'MarkerSize', 12, 'LineWidth', 1, 'color', dat.crossColors{p});
end

bregma = (R * dat.bregma' + t)';
plot(bregma(1), bregma(2), 'x', ...
  'MarkerSize', 12, 'LineWidth', 1, 'color', dat.bregmaColor);
plot(bregma(1), bregma(2), 'o', ...
  'MarkerSize', 8, 'LineWidth', 1, 'color', dat.bregmaColor);
  


transParams.refPts = dat.refPts;
transParams.scaledRefs = scaledRefs;
transParams.selectedPts = dat.selectedPts';
transParams.transSelected = selPts;
transParams.selectedBregma = dat.bregma;
transParams.transBregma = bregma;
transParams.meanRefErr = mean(sqrt(sum((scaledRefs - selPts) .^ 2, 2)));

title(sprintf('Mean ref error: %0.2f', transParams.meanRefErr));


%% Second panel: Aligned image with borders

subplot(1, 3, 2);
showAllenEdgesOnWF(dat.WFIm, transParams, dat.aMaps, 0, 0);

if dat.chkScale.Value  
  title(sprintf('Brain scaled to %0.1f%%', 100 * transParams.scaleConst));
else
  title('No brain scaling allowed');
end


%% Third panel: Aligned image, masked

subplot(1, 3, 3);
showAllenEdgesOnWF(dat.WFIm, transParams, dat.aMaps, 1, 0);

if transParams.meanRefErr < 4
  title('Very nice, well done old chap');
end


%% Button to save transform to workspace

hbw = uicontrol('Style', 'pushbutton', ...
  'String', 'Save params to workspace', ...
  'Units', 'Normalized', 'Position', [0.25 0.05 0.2 0.08], ...
  'ForegroundColor', 'w', ...
  'BackgroundColor', [0 0 1]);
hbw.Callback = @(~, ~) transParamsToWS(transParams, hbw);


%% Button to save transform to opts2.mat, if available

if dat.havePath
  hbs = uicontrol('Style', 'pushbutton', ...
  'String', 'Save params to opts2.mat', ...
  'Units', 'Normalized', 'Position', [0.55 0.05 0.2 0.08], ...
  'ForegroundColor', 'w', ...
  'BackgroundColor', [1 0 0]);
  hbs.Callback = @(~, ~) saveOpts(dat.path, transParams, hbs);
end



function moveCross(src, ev)
% Update data structure after dragging a reference point

hf = gcbf;
dat = hf.UserData;
refI = find(dat.WFCrossHs == src, 1);

% Paranoia, should never happen
if isempty(refI)
  error('Lost track of handle for dragged reference point');
end

dat.selectedPts(refI, :) = [src.XData src.YData];

hf.UserData = dat;


function moveBregma(src, ev)
% Update data structure after dragging a reference point

hf = gcbf;
dat = hf.UserData;

dat.bregma = [src.XData src.YData];

hf.UserData = dat;


function transParamsToWS(transParams, src)

assignin('base', 'transParams', transParams);
src.Enable = 'off';
src.BackgroundColor = [0.5 0.5 1];


function saveOpts(thePath, transParams, src)

loadVar = load(fullfile(thePath, 'opts.mat'));

fnames = fieldnames(loadVar);
if ~isfield(loadVar, 'opts')
  error('opts.mat got deleted!');
end
if length(fnames) > 1
  error('opts.mat contains variables besides opts!');
end

opts = loadVar.opts;
opts.transParams = transParams;

save(fullfile(thePath, 'opts2.mat'), 'opts');

src.Enable = 'off';
src.BackgroundColor = [1 0.5 0.5];
