cPath = 'Y:\data\BpodImager\Animals\';
% cPath = 'X:\smusall\BpodImager\Animals\';
% cPath = 'U:\smusall\BpodImager\Animals\';

iMod = 3;

if iMod == 1
    % GFP data
    [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
    dataOverview = dataOverview(ismember(dataOverview(:,4), 'GFP'), :);
elseif iMod == 2
    % camK data
    [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
    dataOverview = dataOverview(ismember(dataOverview(:,4), 'gcamp'), :);
    cPath = 'U:\smusall\BpodImager\Animals\'; %Widefield data path on grid server
elseif iMod == 3
    % ai93 data
    [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings;
end

% dataOverview = delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
segIdx = {1:54 55:83 84:102 118:135 136:166};

%get allen maps
load('allenDorsalMapSM.mat')
mask = ~dorsalMaps.cutMaskScaled; %mask that is used for all datasets
[x1, y1] = size(mask);
load([cPath animals{1} filesep 'SpatialDisc' filesep recs{1} filesep 'snapshot_1.mat'])
load([cPath animals{1} filesep 'SpatialDisc' filesep recs{1} filesep 'opts2.mat'])
snap = alignAllenTransIm(single(snap),opts.transParams);
[x2, y2] = size(snap);
mask = mask(1:min([x1 x2]), 1:min([y1 y2])); %cut mask to size
rightHs = find(ismember(dorsalMaps.sidesSplit,'L')); %index for labels on the left HSf
allenMask = mask;
halfMask = allenMask(:,1:size(allenMask,2)/2);

%% get unique explained variance
visMap = NaN(sum(~allenMask(:)), length(animals));
visSegMap = NaN(sum(~allenMask(:)), length(animals));
visSegMovie = NaN(sum(~allenMask(:)), 6, length(animals));                
visMovie = NaN(sum(~allenMask(:)), 189, length(animals));                
audMap = NaN(sum(~allenMask(:)), length(animals));
audSegMovie = NaN(sum(~allenMask(:)), 6, length(animals));                
audMovie = NaN(sum(~allenMask(:)), 189, length(animals));                
fullMap = NaN(sum(~allenMask(:)), length(animals));
fullSegMovie = NaN(sum(~allenMask(:)), 6, length(animals));                

for iAnimals = 1 : length(animals)
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep];
    disp(fPath); %current data path
        
    load([fPath 'mask.mat'],'mask');
    load([fPath 'opts2.mat']);
    
    load([fPath 'predVariance' filesep 'newShCurrent' filesep 'fullcorr.mat'],'cMap', 'segMovie', 'cMovie');
    cMap = alignAllenTransIm(arrayShrink(cMap.^2,mask,'split'),opts.transParams);
    segMovie = alignAllenTransIm(arrayShrink(segMovie.^2,mask,'split'),opts.transParams);
    cMovie = alignAllenTransIm(arrayShrink(cMovie.^2,mask,'split'),opts.transParams);
    fullMap(:,iAnimals) = arrayShrink(cMap(:,1:y1),allenMask);
    fullSegMovie(:,:,iAnimals) = arrayShrink(segMovie(:,1:y1,:),allenMask);
    fullMovie = arrayShrink(cMovie(:,1:y1,:),allenMask);
    
    load([fPath 'predVariance' filesep 'newShCurrent' filesep 'visChoicecorr.mat'],'cMap', 'segMovie', 'cMovie');
    cMap = alignAllenTransIm(arrayShrink(cMap.^2,mask,'split'),opts.transParams);
    segMovie = alignAllenTransIm(arrayShrink(segMovie.^2,mask,'split'),opts.transParams);
    cMovie = alignAllenTransIm(arrayShrink(cMovie.^2,mask,'split'),opts.transParams);
    visMap(:,iAnimals) = fullMap(:,iAnimals) - arrayShrink(cMap(:,1:y1),allenMask);
    visSegMovie(:,:,iAnimals) = fullSegMovie(:,:,iAnimals) - arrayShrink(segMovie(:,1:y1,:),allenMask);
    visMovie(:,:,iAnimals) = fullMovie - arrayShrink(cMovie(:,1:y1,:),allenMask);
    
    load([fPath 'predVariance' filesep 'newShCurrent' filesep 'audChoicecorr.mat'],'cMap', 'segMovie', 'cMovie');
    cMap = alignAllenTransIm(arrayShrink(cMap.^2,mask,'split'),opts.transParams);
    segMovie = alignAllenTransIm(arrayShrink(segMovie.^2,mask,'split'),opts.transParams);
    cMovie = alignAllenTransIm(arrayShrink(cMovie.^2,mask,'split'),opts.transParams);
    audMap(:,iAnimals) = fullMap(:,iAnimals) - arrayShrink(cMap(:,1:y1),allenMask);
    audSegMovie(:,:,iAnimals) = fullSegMovie(:,:,iAnimals) - arrayShrink(segMovie(:,1:y1,:),allenMask);
    audMovie(:,:,iAnimals) = fullMovie - arrayShrink(cMovie(:,1:y1,:),allenMask);

end

%%
expChoice = arrayShrink(nanmean(cat(3, visMovie(:,:,visExp), audMovie(:,:,~visExp)),3), allenMask, 'split');
novChoice = arrayShrink(nanmean(cat(3, visMovie(:,:,~visExp), audMovie(:,:,visExp)),3), allenMask, 'split');

expSegChoice = arrayShrink(nanmean(cat(3, visSegMovie(:,:,visExp), audSegMovie(:,:,~visExp)),3), allenMask, 'split');
novSegChoice = arrayShrink(nanmean(cat(3, visSegMovie(:,:,~visExp), audSegMovie(:,:,visExp)),3), allenMask, 'split');

expChoice = arrayShrink(nanmean(cat(2, visMap(:,visExp), audMap(:,~visExp)),2), allenMask, 'split');
novChoice = arrayShrink(nanmean(cat(2, visMap(:,~visExp), audMap(:,visExp)),2), allenMask, 'split');


expChoice = smoothCol(expChoice,10,'box',[],3);
novChoice = smoothCol(novChoice,10,'box',[],3);

aExpBeta = arrayShrink(nanmean(audMovie(:,:,~visExp),3), allenMask, 'split');
aNovBeta = arrayShrink(nanmean(audMovie(:,:,visExp),3), allenMask, 'split');


expSegChoice = arrayShrink(nanmean(cat(3, visSegMovie(:,:,visExp), audSegMovie(:,:,~visExp)),3), allenMask, 'split');
novSegChoice = arrayShrink(nanmean(cat(3, visSegMovie(:,:,~visExp), audSegMovie(:,:,visExp)),3), allenMask, 'split');
clear visMovie audMovie

%% get model weights
visMovie = NaN(sum(~halfMask(:)), 189, length(animals));                
% audMovie = NaN(sum(~halfMask(:)), 189, length(animals));                
for iAnimals = 1 : length(animals)
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep];
    disp(fPath); %current data path
        
    load([fPath 'mask.mat'],'mask');
    load([fPath 'Vc.mat'], 'U');
    load([fPath 'opts2.mat']);
    U = alignAllenTransIm(U,opts.transParams);
    U = arrayShrink(U(:,1:y1,:), allenMask);

    load([fPath 'dimBeta.mat']);
    load([fPath 'regData.mat'], 'idx' ,'trialIdx', 'recIdx', 'recLabels');
    
    cIdx = ismember(recIdx(~idx), find(ismember(recLabels,'Choice')));
    cMovie = U * dimBeta(cIdx, :)';
    cMovie = smoothCol(cMovie,5,'gauss',[],2);
    cMovie = arrayShrink(cMovie,allenMask,'split');
    cMovie = (cMovie(:,1:586/2,:) - cMovie(:,end : -1 : 586/2 + 1,:));
    visMovie(:,1:sum(cIdx),iAnimals) = arrayShrink(cMovie,halfMask);

    
%     cIdx = ismember(recIdx(~idx), find(ismember(recLabels,'audChoice')));
%     cMovie = U * dimBeta(cIdx, :)';
%     cMovie = smoothCol(cMovie,5,'gauss',[],2);
%     cMovie = arrayShrink(cMovie,allenMask,'split');
%     cMovie = (cMovie(:,1:586/2,:) - cMovie(:,end : -1 : 586/2 + 1,:)); 
%     audMovie(:,1:sum(cIdx),iAnimals) = arrayShrink(cMovie,halfMask);
end

%%
% visMovie = smoothCol(visMovie,15,'gauss',[],2);
% audMovie = smoothCol(audMovie,15,'gauss',[],2);

%%
expBeta = arrayShrink(nanmean(cat(3, visMovie(:,:,visExp), audMovie(:,:,~visExp)),3), halfMask, 'split');
novBeta = arrayShrink(nanmean(cat(3, visMovie(:,:,~visExp), audMovie(:,:,visExp)),3), halfMask, 'split');
temp = ((expBeta)*-1) - ((novBeta)*-1);

vExpBeta = arrayShrink(nanmean(visMovie(:,:,visExp),3), halfMask, 'split');
vNovBeta = arrayShrink(nanmean(visMovie(:,:,~visExp),3), halfMask, 'split');

aExpBeta = arrayShrink(nanmean(audMovie(:,:,~visExp),3), halfMask, 'split');
aNovBeta = arrayShrink(nanmean(audMovie(:,:,visExp),3), halfMask, 'split');


temp1 = abs(vExpBeta - vNovBeta);
temp2 = aExpBeta - aNovBeta;

% temp3 = cat(4,(nanmean(cat(4,vExpBeta,aExpBeta),4)),(nanmean(cat(4,vNovBeta,aNovBeta),4)));

%%
% temp1(:,1:586/2,:,:) = (expBeta(:,1:586/2,:,:) - expBeta(:,end : -1 : 586/2 + 1,:,:));
% temp2(:,1:586/2,:,:) = (novBeta(:,1:586/2,:,:) - novBeta(:,end : -1 : 586/2 + 1,:,:));
% temp = abs(temp1) + abs(temp2);

figure;
for iMod = 1:3
    subplot(1,3,iMod);
    temp = nanmean(vExpBeta(:, :, round((segIdxRealign{iMod + 2}))), 3)*-1;
%     temp = nanmean(expBeta(:, :, round((segIdx{iMod + 2}))) * -1, 3) - nanmean(novBeta(:, :, round((segIdx{iMod + 2}))) * -1, 3);
    mapImg = imagesc(smooth2a(temp, 5, 5)); axis image;
    hold on; caxis([-0.0015 0.0015]);
    for x = 1: length(rightHs)
        plot(dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,2), dorsalMaps.edgeOutlineSplit{rightHs(x)}(:,1), 'w', 'LineWidth', 0.2);axis image
    end
    set(mapImg,'AlphaData',~isnan(mapImg.CData));
end
colormap inferno