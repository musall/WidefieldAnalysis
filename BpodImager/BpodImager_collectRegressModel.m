function allData = BpodImager_collectRegressModel(cMod,dType,reload,verbose)

%% select data sets
dataOverview = delayDecRecordings; %get overview for all recordings

if ~exist('reload','var') || isempty(reload)
    reload = false;  % Don't reload data from server of not required
end

if ~exist('verbose','var') || isempty(verbose)
    verbose = false;  % flag to show additional feedback. Shows which data is used and plots vessel images.
end

%% general variables
Paradigm = 'SpatialDisc';
cPath = 'X:\smusall\BpodImager\Animals\'; %Widefield data path on grid server
sPath = 'U:\space_managed_data\BpodImager\Animals\'; %Widefield data path on grid server
pxPerMM = 206/4;    % pixels per milimeters - default is 206/4
rotAngle = 40;      % Angle to rotate imaging data

%%
animals = dataOverview(:,1);
Cnt = 0;

for iAnimals = 1:length(animals)
    if strcmpi(dataOverview{iAnimals,2},cMod) || strcmpi(cMod,'all')
        Cnt = Cnt +1;

        %% load data
        fPath = [cPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
        if isempty(dir([fPath 'fullBeta.mat'])) || reload %check if file exists on hdd and pull from server otherwise
            if ~isdir(fPath); mkdir(fPath); end
            sourcePath = [sPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
            copyfile([sourcePath 'fullBeta.mat'],[fPath 'fullBeta.mat']);
            copyfile([sourcePath 'mask.mat'],[fPath 'mask.mat']);
            copyfile([sourcePath 'opts2.mat'],[fPath 'opts2.mat']);
        end
        load([fPath 'Snapshot_1.mat']);
        load([fPath 'mask.mat'])
        load([fPath 'opts2.mat'])
        
        if isempty(dir([fPath dType 'B.mat'])) || reload %check if requested regressor set is found
            load([fPath 'fullBeta.mat'])
            idx = recIdx == find(ismember(recLabels,dType));
            
            if sum(idx) == 0
                disp('Requested beta set not found. Try using one of the following labels.')
                disp('===============================')
                disp(recLabels')
                disp('===============================')
                return
            else
                data = fullBeta(:,:,idx);
                save([fPath dType 'B.mat'], 'data', 'recLabels','-v7.3'); %save current beta map for next time
                clear fullBeta
            end
        else
            load([fPath dType 'B.mat']); %load current beta map
            data = single(reshape(data,size(mask,1),size(mask,2),[]));
        end
        
        %% collect data
        if Cnt == 1
            % create mask from allen coordinates
            load('allenDorsalMap.mat')
            addEdgeOutlinesToDorsalMap(dorsalMaps); %make area figure for orientation
            allenMask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
            [x1, y1] = size(allenMask);
            [x2, y2] = size(data);
            allenMask = allenMask(1:min([x1 x2]), 1:min([y1 y2])); %cut allen mask to size
            [x1, y1] = size(allenMask);
        end
        
        data = alignAllenTransIm(data,opts.transParams);
        snap = alignAllenTransIm(snap,opts.transParams);
        combData{Cnt} = arrayShrink(arrayShrink(data(1:x1,1:y1,:),allenMask),allenMask,'split');
        combSnap{Cnt} = arrayShrink(arrayShrink(snap(1:x1,1:y1),allenMask),allenMask,'split');
        
    end
end

%% show vessel images for selected animals and get largest image dimensions
if strcmpi(cMod,'all')
    idx = 1:size(dataOverview,1);
else
    idx = strcmpi(dataOverview(:,2),cMod);
end
animals = dataOverview(idx, 1);
recs = dataOverview(idx, 3);
maxX = 0;
maxY = 0;

for x = 1:Cnt
    if verbose
        if x == 1; figure; end
        subplot(2,ceil(Cnt/2),x);
        imagesc(combSnap{x});axis image; colormap gray;
        hline(size(combSnap{x},1)/2,'--w')
        vline(size(combSnap{x},2)/2,'--w')
        title([animals{x} ' - ' recs{x}])
    end
    temp = size(combData{x});
    if temp(1) > maxY
        maxY = temp(1);
    end
    if temp(2) > maxX
        maxX = temp(2);
    end
end

%% combine aligned individual data
allData = NaN(maxY ,maxX , size(combData{1},3), Cnt, 'single');

for x = 1:Cnt
    temp = size(combData{x});
    shift(1) = (maxY - temp(1)) / 2;
    shift(2) = (maxX - temp(2)) / 2;
    allData(shift(1) + (1 : temp(1)), shift(2) + (1 : temp(2)) ,:, x) = single(combData{x});
end

%% give some feedback
labels = dataOverview(strcmpi(dataOverview(:,2),cMod),1);
for x = 1:length(labels)
    labels{x} = [labels{x} ' - ' dType];
end

if nargout == 0
    compareMovie(allData, labels);
    allData = [];
end