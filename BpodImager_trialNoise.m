% function BpodImager_trialNoise

dataOverview = delayDecRecordings;
animals = dataOverview(:,1)';
Paradigm = 'SpatialDisc';
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
timeVar = NaN(sum(~allenMask(:)),3,length(animals),'single');
trialVar = NaN(sum(~allenMask(:)),3,length(animals),'single');

%%
for iAnimals = 1 : size(dataOverview,1)

    fPath = [cPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
    load([fPath 'Vc.mat'], 'U'); %load spatial data
    load([fPath 'mask.mat'], 'mask'); %load mask
    load([fPath 'opts2.mat'], 'opts'); %load alignment
    load([fPath 'interpVc.mat'], 'Vc', 'frames'); %load real data
    load([fPath 'interpVtask.mat'], 'Vtask'); %load model reconstruction data
    load([fPath 'interpVspontMotor.mat'], 'VspontMotor'); %load model reconstruction data
    load([fPath 'interpVopMotor.mat'], 'VopMotor'); %load model reconstruction data
    U = arrayShrink(U,mask,'merge');
    
    %reshape to trialized format
    Vc = reshape(Vc,size(Vc,1),frames,[]);
    Vtask = reshape(Vtask,size(Vtask,1),frames,[]);
    VopMotor = reshape(VopMotor,size(VopMotor,1),frames,[]);
    VspontMotor = reshape(VspontMotor,size(VspontMotor,1),frames,[]);
    
    % compute trial-by-trial and trial-averaged correlations
    [trialVar(:,1,iAnimals) , timeVar(:,1,iAnimals)] = rateDisc_computeTrialCorr(U, Vc, Vtask, mask, allenMask, opts);
    [trialVar(:,2,iAnimals) , timeVar(:,2,iAnimals)] = rateDisc_computeTrialCorr(U, Vc, VopMotor, mask, allenMask, opts);
    [trialVar(:,3,iAnimals) , timeVar(:,3,iAnimals)] = rateDisc_computeTrialCorr(U, Vc, VspontMotor, mask, allenMask, opts);
    
end

%% make some figures
titles = {'Task' 'Instructed' 'Uninstructed'};
figure;
temp = arrayShrink(nanmean(trialVar,3),allenMask,'split');
meanData = squeeze(mean(nanmean(trialVar,1),3));
semData = squeeze(sem(nanmean(trialVar,1),3));
cRange = [0 0.4];

for x = 1 : 3
    subplot(1,4,x);
    mapImg = imshow(temp(:,:,x),cRange);
    colormap(mapImg.Parent,'inferno'); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['TrialByTrial - ' titles{x}]); colorbar;
end
subplot(1,4,4);
errorbar(meanData,semData,'ok','linewidth',2,'MarkerFaceColor','w','MarkerSize',8);
axis square; xlim([0.5 length(meanData)+0.5]);
title('Pixel average'); ylim([0 0.4]);

figure;
temp = arrayShrink(nanmean(timeVar,3),allenMask,'split');
meanData = squeeze(mean(nanmean(timeVar,1),3));
semData = squeeze(sem(nanmean(timeVar,1),3));
cRange = [0.5 1];

for x = 1 : 3
    subplot(1,4,x);
    mapImg = imshow(temp(:,:,x),cRange);
    colormap(mapImg.Parent,'inferno'); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['TrialAverage - ' titles{x}]); colorbar;
end
subplot(1,4,4);
errorbar(meanData,semData,'ok','linewidth',2,'MarkerFaceColor','w','MarkerSize',8);
axis square; xlim([0.5 length(meanData)+0.5]);
title('Pixel average'); ylim([0.8 1]);




