function [allData, combInd, labels] = BpodImager_compareModalities(cMod,dType,regType,noPlot)

if ~exist('regType', 'var')
    regType = [];
end

if ~exist('noPlot', 'var')
    noPlot = false;
end

if ~isempty(regType)
    regType = ['-' regType 'Rebuild'];
end

%% basic variables
Paradigm = 'SpatialDisc';
cPath = 'H:\BpodImager\Animals\'; %Widefield data path
sPath = 'U:\space_managed_data\\BpodImager\Animals\'; %Widefield data path on server
trialSegments = [60 90 141 171];
segLabels = {'Lever' 'Stimulus' 'Wait' 'Response'};

if contains(dType,'Expert')
    trainMod = true; %flag to invert dPrime sign to compare dPrime for trained modality instead of vision-audio
    dType = strrep(dType,'Expert','');
else
    trainMod = false;
end

%%
dataOverview = delayDecRecordings; %get overview for all recordings
animals = dataOverview(:,1);
Cnt = 0;

for iAnimals = 1:length(animals)
    if strcmpi(dataOverview{iAnimals,2},cMod) || strcmpi(cMod,'All')
        
        fPath = [cPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
        load([fPath animals{iAnimals} '_opts.mat']);
        opts.modality = dataOverview{iAnimals,2};
        save([fPath animals{iAnimals} '_opts.mat'],'opts');
        
        Cnt = Cnt+1;
        allPxPerMM(Cnt) = opts.pxPerMM; %keep pxPerMM to compare across datasets
        disp(fPath)
        load([fPath animals{iAnimals} '-' dataOverview{iAnimals,3} '-psthVision'],'newMask');
        load([fPath 'Snapshot_1.mat']);
        snap = Widefield_mapAlign(snap,opts); %rotate for straight midline and set bregma to center
        
        if contains(dType,'Prime') %if asking for dPrime data
            data = load([fPath filesep animals{iAnimals} '-' dataOverview{iAnimals,3} regType '-AVmodPrime.mat'],dType);
        elseif contains(dType,'aud') %if asking for auditory PSTH
            data = load([fPath filesep animals{iAnimals} '-' dataOverview{iAnimals,3} regType '-psthAudio.mat'],dType);
        elseif contains(dType,'vis') %if asking for visual PSTH
            data = load([fPath filesep animals{iAnimals} '-' dataOverview{iAnimals,3} regType '-psthVision.mat'],dType);
        elseif contains(dType,'corr','IgnoreCase',true) && ~contains(dType,'Prime') %if asking for predicted variance data, make sure corresponding file is present
            if exist([fPath dType '.mat'],'file') ~= 2
                copyfile([strrep(fPath,cPath,sPath) dType '.mat'],[fPath dType '.mat']);
                copyfile([strrep(fPath,cPath,sPath) 'mask.mat'],[fPath 'mask.mat']);
            end
            data = load([fPath dType '.mat'],'cMap');
            load([fPath 'mask.mat'],'mask');
            data.cMap = Widefield_mapAlign(arrayShrink(data.cMap,mask,'split'),opts);
            data.cMap = arrayShrink(data.cMap,newMask);
        end
        
        if isempty(fieldnames(data))
            error('Error. Requested data type not found');
        else
            cField = fieldnames(data);
            singleData = data.(cField{1});
            if trainMod && contains(dType,'Prime') && strcmpi(dataOverview{iAnimals,2},'Audio')
                singleData = -singleData; %invert dPrime sign for audio experts
            end
            singleData = arrayShrink(singleData,newMask,'split');
            
            if contains(dType,'sel') %if asking for selected trials
                indData = load([fPath filesep animals{iAnimals} '-' dataOverview{iAnimals,3} '-trialIdx.mat'],'selInd');
            elseif contains(dType,'vis') %if asking for visual trials
                indData = load([fPath filesep animals{iAnimals} '-' dataOverview{iAnimals,3} '-trialIdx.mat'],'vInd');
            elseif contains(dType,'aud') %if asking for audio trials
                indData = load([fPath filesep animals{iAnimals} '-' dataOverview{iAnimals,3} '-trialIdx.mat'],'aInd');
            elseif contains(dType,'corrPrime') %if asking for succesful trials
                indData = load([fPath filesep animals{iAnimals} '-' dataOverview{iAnimals,3} '-trialIdx.mat'],'sucInd');
            elseif contains(dType,'allPrime') || strcmpi(cField,'cMap') %if asking for all trials
                indData = load([fPath filesep animals{iAnimals} '-' dataOverview{iAnimals,3} '-trialIdx.mat'],'allInd');
            end
            cField = fieldnames(indData);
            indData = indData.(cField{1});
        end
        
        allSnap{Cnt} = snap;
        allMask{Cnt} = newMask;
        combData{Cnt} = singleData;
        combInd{Cnt} = indData;
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

if ~noPlot
    figure
end
    
for x = 1:Cnt
    
    if ~noPlot
        subplot(2,ceil(Cnt/2),x);
        imagesc(allSnap{x});axis image; colormap gray;
        hline(size(allSnap{x},1)/2,'--w')
        vline(size(allSnap{x},2)/2,'--w')
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

%% combine data
allData = NaN(maxY ,maxX , size(combData{1},3), Cnt, 'single');
for x = 1:Cnt
    temp = size(combData{x});
    shift(1) = (maxY - temp(1)) / 2;
    shift(2) = (maxX - temp(2)) / 2;
    allData(shift(1) + (1 : temp(1)), shift(2) + (1 : temp(2)) ,:, x) = single(combData{x});
end

% if size(allData,3) > 1
%     dataAvg = nanmean(allData(:,:,1:30,:),3);
%     allData = bsxfun(@minus, allData, dataAvg); % correct baseline offset (usually Vc has zero mean)
% end

%% show averages for trial segments
if ~noPlot
    figure('name','Average maps');
    cData = nanmean(allData,4); %combine individual animals
    dataAvg = nanmean(cData(:,:,1:30),3);
    cData = bsxfun(@minus, cData, dataAvg); % correct baseline offset
    for iSegs = 1:length(trialSegments)
        
        idx = trialSegments(iSegs) : trialSegments(iSegs) + 15;
        subplot(2,ceil(length(trialSegments)/2),iSegs);
        cImg = imshow(nanmean(cData(:,:,idx),3),[-nanstd(cData(:)) nanstd(cData(:))]*2); cmap = colormap_blueblackred;
        colormap(cImg.Parent,cmap)
        axis image; colorbar;
        set(cImg,'AlphaData',~isnan(cData(:,:,1))); %make NaNs transparent.
        title(segLabels{iSegs})
        
    end
end
%% show full movies
if strcmpi(cMod,'all')
    labels = dataOverview(1:size(allData,4),1);
else
    labels = dataOverview(strcmpi(dataOverview(:,2),cMod),1);
end
for iLabels = 1:size(allData,4)
    labels{iLabels} = [labels{iLabels} ' - ' dType regType ' - ' num2str(sum(combInd{iLabels})) '\' num2str(length(combInd{iLabels}))];
end

if nargout == 0
    compareMovie(allData, labels);
    allData = [];
end