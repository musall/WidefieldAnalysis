% function BpodImager_peerPredict

pxMerge = 1; %merge 'pxMerge' pixels into blocks

dataOverview = delayDecRecordings;
animals = dataOverview(:,1)';
Paradigm = 'SpatialDisc';
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
load('allenDorsalMapSM.mat')
allenMask = dorsalMaps.allenMask;
[xRange, yRange] = rateDisc_maskRange(allenMask,500); % get inner range of allenMask
allenMask = allenMask(yRange,xRange);
allenMask = arrayResize(allenMask,pxMerge) == 1;

%%
for iAnimals = 1 : size(dataOverview,1)

    fPath = [cPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
    load([fPath 'Vc.mat'], 'U'); %load spatial data
    load([fPath 'mask.mat'], 'mask'); %load mask
    load([fPath 'opts2.mat'], 'opts'); %load alignment
    load([fPath 'interpVc.mat'], 'Vc', 'frames'); %load real data
    
    rng default %reset randum number generator
    randIdx = randperm(size(Vc,2)); %generate randum number index if required
    fitIdx = randIdx(1:round(size(Vc,2)/2));
    testIdx = randIdx(round(size(Vc,2)/2)+1:end);
    
    % resize and isolate pixels
    U = arrayResize(U,pxMerge);
    mask = arrayResize(mask,pxMerge) ~= 0;
    U = arrayShrink(U,mask);
    
    % compute block-dataset and predict blocks
    pxRand = randperm(size(U,1));
    cData = U(pxRand(1:500),:)*Vc;
    for x = 1:size(cData,1)
        cIdx = true(size(cData,1),1); cIdx(x) = false;
        [ridgeVals, dimBeta] = ridgeMML(cData(~cIdx,fitIdx)', cData(cIdx,fitIdx)', true, 1); %get ridge penalties and beta weights.
        
        a = cData(~cIdx,testIdx)';
        b = cData(cIdx,testIdx)' * dimBeta;
        allFit(x) = corr2(a,b);
    end
    
    

end

%% make some figures
figure;
ax = plot(-frames:frames,autoSpec, 'color',[0.5 0.5 0.5]); axis square; hold on;
plot(-frames:frames,mean(autoSpec,2), 'k', 'linewidth',4)
xlim([-frames frames]);
xlabel('Frame shift'); ylabel('Correlation coefficient');

cData = double(mean(autoSpec(frames+1:frames+100,:),2));
% cData = exp(-(1:100)/10)';
x = double(1:length(cData))';
g = fittype('a-b*exp(-x/c)');
f0 = fit(x,cData,g,'StartPoint',[[ones(size(x)), -exp(-x)]\cData; 1]);
title(['Vc autocorrelation spectrum; Tau = ' num2str((f0.c * (1/30))*1000) 'ms']);
vline([-f0.c f0.c]);

%% autocorrelation / model fit maps
figure;
temp1 = arrayShrink(nanmean(autoVar,3),allenMask,'split');
temp2 = arrayShrink(nanmean(modelFit,3),allenMask,'split');
cRange = [0 0.8];
for x = 1 : size(autoVar,2)
    subplot(2,size(autoVar,2),x);
    mapImg = imshow(temp1(:,:,x),[0 0.7]);
    colormap(mapImg.Parent,'inferno'); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['R^2, AutoCorr - ' num2str(30/freqChange(x)) 'Hz']); colorbar;
    
    subplot(2,size(autoVar,2),size(autoVar,2) + x);
    mapImg = imshow(temp2(:,:,x),[0 0.7]);
    colormap(mapImg.Parent,'inferno'); axis image
    set(mapImg,'AlphaData',~isnan(mapImg.CData)); %make NaNs transparent.
    title(['R^2, AutoCorr - ' num2str(30/freqChange(x)) 'Hz']); colorbar;
end

%% errorbar plot
clear meanData semData
meanData(1,:) = squeeze(mean(nanmean(autoVar,1),3));
meanData(2,:) = squeeze(mean(nanmean(modelFit,1),3));
semData(1,:) = squeeze(sem(nanmean(autoVar,1),3));
semData(2,:) = squeeze(sem(nanmean(modelFit,1),3));

figure;
errorbar((1:3)-0.15, meanData(1,:),semData(1,:),'.k','linewidth',2); hold on;
errorbar((1:3)+0.15, meanData(2,:),semData(2,:),'.k','linewidth',2);
ax = bar(meanData'); ax = ax.Parent; hold off
axis square
set(gca,'xTick',1:length(meanData))
set(gca,'xTickLabel',num2cell(30./freqChange))
ylim([0 0.75]); xlabel('Sampling rate (Hz)'); ylabel('Correlation');
legend(ax.Children(2:-1:1),{'AutoCorrelation' 'Full model'});
title('Relation between full model and autocorrelation')
