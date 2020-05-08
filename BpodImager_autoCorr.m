% function BpodImager_autoCorr

dataOverview = delayDecRecordings;
animals = dataOverview(:,1)';
Paradigm = 'SpatialDisc';
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
gPath = '\\grid-hs\\churchland_hpc_home\smusall\BpodImager\Animals\'; %data path on the server
load('allenDorsalMapSM.mat')

allenMask = dorsalMaps.allenMask;
autoVar = NaN(sum(~allenMask(:)),3,length(animals),'single');
modelFit = NaN(sum(~allenMask(:)),3,length(animals),'single');

autoShift = 45; % number of pixels for shift
autoSpec = NaN(autoShift*2+1,length(animals),'single');
modelShift = NaN(autoShift*2+1,length(animals),'single');
crossSpec = NaN(autoShift*2+1,length(animals),'single');
freqChange = [1 10 30]; %downsample frequencies to 30, 3, and 1Hz

%%
for iAnimals = 1 : size(dataOverview,1)

    fPath = [cPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
    gfPath = [gPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep 'shiftVals' filesep];
    load([fPath 'Vc.mat'], 'U'); %load spatial data
    load([fPath 'mask.mat'], 'mask'); %load mask
    load([fPath 'opts2.mat'], 'opts'); %load alignment
    load([fPath 'interpVc.mat'], 'Vc', 'frames'); %load real data
    load([fPath 'interpVfull.mat'], 'Vm'); %load predicted data
    U = arrayShrink(U,mask,'merge');
    
    Vc = reshape(Vc,size(Vc,1),frames,[]);
    Vm = reshape(Vm,size(Vm,1),frames,[]);
    Vc(:,[1 end],:) = 0; %set first and last frame to 0
    Vm(:,[1 end],:) = 0; %set first and last frame to 0
    Vc = reshape(Vc,size(Vc,1),[]);
    Vm = reshape(Vm,size(Vm,1),[]);
    
    % compute autocorrelation spectrum
    a = zeros(size(Vc,1),autoShift*2+1);
    b = zeros(size(Vc,1),autoShift*2+1);
    for x = 1 : size(Vc,1)
        a(x,:) = xcorr(Vc(x,:),Vc(x,:),autoShift);
        b(x,:) = xcorr(Vc(x,:),Vm(x,:),autoShift);
    end
    autoSpec(:,iAnimals) = mean(a)./max(mean(a));
    crossSpec(:,iAnimals) = mean(b)./max(mean(b));
    
    % load model shift results
    Cnt = 0;
    for iShifts = autoShift : -1 : -autoShift %do this in reverse because neural data was shifted against behavior here.
        Cnt = Cnt + 1;
        load([gfPath 'shiftCorr_' num2str(iShifts)])
        modelShift(Cnt, iAnimals) = corrMat;
        xVals(Cnt) = -iShifts;
    end
    
    for x = 1 : length(freqChange)
        % resample data to 15Hz
        lowVc = single(resample(double(Vc)', 1, freqChange(x))');
        lowVm = single(resample(double(Vm)', 1, freqChange(x))');
        
        % compute auto-correlations after shifting by 1 frame. This is to identify the non-white noise part of the data
        corrMat = rateDisc_modelCorr(lowVc,[zeros(size(lowVc,1),1) lowVc(:,1:end-1)],U) .^2; %compute explained variance after shifting data by one frame
        corrMat = arrayShrink(corrMat,mask,'split');
        corrMat = alignAllenTransIm(corrMat,opts.transParams);
        autoVar(:,x,iAnimals) = arrayShrink(corrMat(:,1:size(allenMask,2)),allenMask, 'merge');
        
        % compute correlation between full model and data again. Just to make sure the numbers come out the same.
        corrMat = rateDisc_modelCorr(lowVc,lowVm,U) .^2; %compute explained variance after shifting data by one frame
        corrMat = arrayShrink(corrMat,mask,'split');
        corrMat = alignAllenTransIm(corrMat,opts.transParams);
        modelFit(:,x,iAnimals) = arrayShrink(corrMat(:,1:size(allenMask,2)),allenMask, 'merge');
    end
    fprintf('Done recording %d/%d \n',iAnimals,size(dataOverview,1))
end

%% make some figures
figure;
subplot(1,3,1);
ax = plot(-autoShift:autoShift,autoSpec.^2, 'color',[0.5 0.5 0.5]); axis square; hold on;
plot(-autoShift:autoShift,mean(autoSpec.^2,2), 'k', 'linewidth',4)
xlim([-autoShift autoShift]);
xlabel('Frame shift'); ylabel('Explaiend variance');

cData = double(mean(autoSpec(autoShift+1:autoShift+autoShift,:).^2,2));
% cData = exp(-(1:100)/10)';
x = double(1:length(cData))';
g = fittype('a-b*exp(-x/c)');
f0 = fit(x,cData,g,'StartPoint',[[ones(size(x)), -exp(-x)]\cData; 1]);
title(['Vc autocorrelation spectrum; Tau = ' num2str((f0.c * (1/30))*1000) 'ms']);
ylim([-0.2 1.2]); vline([-f0.c f0.c]);
xlim([-autoShift autoShift]);

subplot(1,3,2);
ax = plot(-autoShift:autoShift,crossSpec.^2, 'color',[0.5 0.5 0.5]); axis square; hold on;
plot(-autoShift:autoShift,mean(crossSpec.^2,2), 'k', 'linewidth',4)
xlim([-autoShift autoShift]);
xlabel('Frame shift'); ylabel('Explaiend variance');

cData = double(mean(crossSpec(autoShift+1:autoShift+autoShift,:).^2,2));
% cData = exp(-(1:100)/10)';
x = double(1:length(cData))';
g = fittype('a-b*exp(-x/c)');
f0 = fit(x,cData,g,'StartPoint',[[ones(size(x)), -exp(-x)]\cData; 1]);
title(['Vc autocorrelation spectrum; Tau = ' num2str((f0.c * (1/30))*1000) 'ms']);
ylim([-0.2 1.2]); vline([-f0.c f0.c]);
xlim([-autoShift autoShift]);

subplot(1,3,3);
ax = plot(xVals,modelShift, 'color',[0.5 0.5 0.5]); axis square; hold on;
plot(-autoShift:autoShift,mean(modelShift,2), 'k', 'linewidth',4)
xlim([-autoShift autoShift]);
xlabel('Behavior shift'); ylabel('cvR^2');

cData = double(mean(modelShift(autoShift+1:autoShift+autoShift,:),2));
% cData = exp(-(1:100)/10)';
x = double(1:length(cData))';
g = fittype('a-b*exp(-x/c)');
f0 = fit(x,cData,g,'StartPoint',[[ones(size(x)), -exp(-x)]\cData; 1]);
title(['Shifted behavioral model; Tau = ' num2str((f0.c * (1/30))*1000) 'ms']);
ylim([0 0.45]); vline([-f0.c f0.c]);
xlim([-autoShift autoShift]);


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


%same thing as boxplot
figure
subplot(1,2,1);
boxplot(squeeze(nanmean(autoVar,1))');
set(gca,'xTickLabel',num2cell(30./freqChange))
ylim([0 1]); xlabel('Sampling rate (Hz)'); ylabel('Correlation');

title('AutoCorrelation');
subplot(1,2,2);
boxplot(squeeze(nanmean(modelFit,1))');
title('Model');
set(gca,'xTickLabel',num2cell(30./freqChange))
ylim([0 1]); xlabel('Sampling rate (Hz)'); ylabel('Correlation');