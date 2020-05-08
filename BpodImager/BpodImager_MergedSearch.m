function BpodImager_MergedSearch(Animal)
% Code to analyze SVD compressed widefield data by combining multiple
% recordings from the same animal.
% Not working yet.

%% basic variables
Paradigm = 'SpatialDisc';
bPath = ['H:\BpodImager\Animals\' Animal filesep Paradigm filesep]; %Widefield data path
sPath = ['U:\space_managed_data\BpodImager\Animals\' Animal filesep Paradigm filesep ]; %server data path
smth = 2;           %2-D smooth for plots/movies
rotAngle = 40;       % Angle to rotate imaging data
dimCnt = 250;       % number of dimensions used for analysis
tBin = 1;           % timebins for analysis. Used to average frames together for speed.
pxPerMM = 206/4;    % pixels per milimeters - default is 206/4

%% check recordings for merge
recs = dir(bPath);
recs = recs([recs.isdir] & ~strncmpi('.', {recs.name}, 1));
disp(['Animal ' Animal '; Found ' int2str(length(recs)) ' recordings'])
disp('=================================')

%% load SVD data
for iRecs = 1:length(recs)
    cPath = [bPath recs(iRecs).name filesep];
    disp(cPath)
    
    if isempty(dir([cPath 'Vc.mat'])) %check if file exists on hdd and pull from server otherwise
        copyfile([sPath recs(iRecs).name filesep 'Vc.mat'],[cPath 'Vc.mat']);
        copyfile([sPath recs(iRecs).name filesep 'mask.mat'],[cPath 'mask.mat']);
    end
    
    load([cPath 'Vc.mat'])
    load([cPath 'mask.mat'])
    
    allV{iRecs} = Vc;
    allU{iRecs} = U;
    allTrials{iRecs} = trials;
    allMask{iRecs} = mask;
    clear Vc U
    
end

%% check opts file, create one if not present
if isempty(dir([cPath '\' Animal '_opts.mat'])) %check if opts file exist
    check = true;
    
    while check
        
        opts.modality = questdlg('Select animal modality','Animal expertise','Visual','Audio','Mixed','Visual');
        opts.pxPerMM = pxPerMM;
        
        copyfile([sPath '\Snapshot_1.mat'],[cPath '\Snapshot_1.mat']);
        load([cPath '\Snapshot_1.mat']);
        snap = double(imrotate(snap,rotAngle,'crop'));
        
        figure(90);
        imagesc(snap);axis image; colormap gray; hold on
        
        [bregma(1), bregma(2)] = Widefield_SelectCoordinates(snap,pxPerMM,'bregma'); %set bregma coordinate
        plot(bregma(1), bregma(2), 'or', 'MarkerFaceColor', 'r')
        snap(round(bregma(2)), round(bregma(1))) = inf;
        
        [lambda(1), lambda(2)] = Widefield_SelectCoordinates(snap,pxPerMM,'lambda'); %set lambda coordinate
        plot(lambda(1), lambda(2), 'ob', 'MarkerFaceColor', 'b')
        snap(round(lambda(2)), round(lambda(1))) = inf;
        
        [frontal(1), frontal(2)] = Widefield_SelectCoordinates(snap,pxPerMM,'frontal'); %set that other coordinate
        plot(frontal(1), frontal(2), 'ow', 'MarkerFaceColor', 'w')
        snap(round(frontal(2)), round(frontal(1))) = inf;
        
        degAngle(1) = (atan2(bregma(2) - lambda(2), bregma(1) - lambda(1))*180/pi) * -1;    %angle for lambda to bregma line, relative to x-axis
        degAngle(2) = (atan2(frontal(2) - bregma(2), frontal(1) - bregma(1))*180/pi) * -1;  %angle for bregma to frontal line, relative to x-axis
        imAngle = 90 - mean(degAngle);
        
        snap = imrotate(snap,imAngle,'crop');
        opts.rotAngle = rotAngle + imAngle; %rotation angle of raw data to straighten image
        
        [b, a] = ind2sub(size(snap),find(isinf(snap))); %find landmark coordinates from rotated image
        [b,yy] = sort(b);
        a = a(yy);
        
        opts.frontal = [a(1) b(1)];
        opts.bregma = [a(2) b(2)];
        opts.lambda = [a(3) b(3)];
        
        figure(90); hold off; 
        imagesc(snap);axis image; colormap gray; hold on
        grid(gca,'on');grid minor;set(gca,'GridColor','w');
        plot(opts.bregma(1), opts.bregma(2), 'ob', 'MarkerFaceColor', 'b');
        plot(opts.lambda(1), opts.lambda(2), 'or', 'MarkerFaceColor', 'r');
        plot(opts.frontal(1), opts.frontal(2), 'ow', 'MarkerFaceColor', 'w')
        title('Aligned image');
        
        happy = questdlg('Happy?','Happy?','Yes','No','No');
        
        if strcmpi(happy,'Yes')
            check = false;
        end
    end
    
    save([cPath '\' Animal '_opts.mat'], 'opts')
    
else
    load([cPath '\' Animal '_opts.mat']);
end

%% load behavior data and get indices for different trial types
cFile = ls([cPath filesep Animal '_' Paradigm '*.mat']);
if isempty(cFile)
    sFile = ls([sPath filesep Animal '_' Paradigm '*.mat']);
    copyfile([sPath filesep sFile],[cPath filesep sFile]);
    cFile = ls([cPath filesep Animal '_' Paradigm '*.mat']);
end
load([cPath filesep cFile]); %load behavior data

ind = trials > SessionData.nTrials; %check for imaging trials that are not in the bhv file
trials(ind) = [];

SummaryFigure(SessionData,[cPath filesep Animal '-' Rec '-BehaviorOverview.pdf']); %behavior summary

%% compute delta F/F movies
trace = svdFrameReconstruct(U(:,:,1:dimCnt), mean(Vc(1:dimCnt,:,:),3));  %average movie over all trials
Widefield_SaveToAvi(trace,[Animal '-' Rec '-AllAvg-fps10'],10,'colormap_blueblackred',[-0.02 0.02],cPath,smth,2,opts.rotAngle);

cMod = {'lVision' 'lAudio' 'lAudioVisual';'rVision' 'rAudio' 'rAudioVisual'};
for iMod = 1:3 % StimType = 1 (Vision), StimType = 2 (Audio), StimType = 3 (Audiovisual),
    for iSide = 1:2
        
        detectInd = SessionData.StimType(trials) == iMod & ~SessionData.DidNotChoose(trials) ...
            & ~SessionData.DidNotLever(trials) & SessionData.CorrectSide(trials) == iSide;
        trace = svdFrameReconstruct(U,median(Vc(:,:,detectInd),3)); %get average data
        trace = bsxfun(@minus,trace,  mean(trace(:,:,1:20),3)); %correct for baseline
        Widefield_SaveToAvi(trace,[Animal '-' Rec '-' cMod{iSide,iMod} '-fps10'],10,'colormap_blueblackred',[-0.01 0.01],cPath,smth,2,opts.rotAngle);

    end
end

%% select trials for AUC
selTrials = Behavior_AnalysisTrialselection(SessionData,trials,1,2,true); %select for even number of singe modality trials and choices for each side
selInd = false(1,length(trials));
selInd(selTrials) = true;

selErrorTrials = Behavior_AnalysisTrialselection(SessionData,trials,1,2,false); %select for even number of singe modality trials and choices for each side
selErrorInd = false(1,length(trials));
selErrorInd(selErrorTrials) = true;

%% get trial indices
vInd = SessionData.StimType(trials) == 1; % visual trials
aInd = SessionData.StimType(trials) == 2; % audio trials
allInd = vInd | aInd; %all unisensory trials
sucInd = SessionData.Rewarded(trials) & SessionData.Assisted(trials) & allInd; %find succesful unisensory trials

save([cPath filesep Animal '-' Rec '-vInd.mat'],'vInd');
save([cPath filesep Animal '-' Rec '-aInd.mat'],'aInd');
save([cPath filesep Animal '-' Rec '-allInd.mat'],'allInd');
save([cPath filesep Animal '-' Rec '-sucInd.mat'],'sucInd');
save([cPath filesep Animal '-' Rec '-selInd.mat'],'selInd');
save([cPath filesep Animal '-' Rec '-selErrorInd.mat'],'selErrorInd');

%% compute averages and save
sucVis = svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,sucInd & vInd),3))); %reconstruct frames
save([cPath filesep Animal '-' Rec '-sucVis.mat'],'sucVis');
clear sucVis

selVis = svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,selInd & vInd),3))); %reconstruct frames
save([cPath filesep Animal '-' Rec '-selVis.mat'],'selVis');
clear selVis

selErrorVis = svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,selErrorInd & vInd),3))); %reconstruct frames
save([cPath filesep Animal '-' Rec '-selErrorVis.mat'],'selErrorVis');
clear selErrorVis

sucAudio = svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,sucInd & aInd),3))); %reconstruct frames
save([cPath filesep Animal '-' Rec '-sucAudio.mat'],'sucAudio');
clear sucAudio

selAudio = svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,selInd & aInd),3))); %reconstruct frames
save([cPath filesep Animal '-' Rec '-selAudio.mat'],'selAudio');
clear selAudio

selErrorAudio = svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,selErrorInd & aInd),3))); %reconstruct frames
save([cPath filesep Animal '-' Rec '-selErrorAudio.mat'],'selErrorAudio');
clear selErrorVis

%% AUC analysis - Modality difference
frameCnt = size(Vc,2);
smallMask = arrayResize(mask,2) == 1;
% time = ((1/sRate)*tBin:(1/sRate)*tBin:4)'-2;

aAucMap = zeros([sum(~smallMask(:)),frameCnt/tBin]);
cAucMap = zeros([sum(~smallMask(:)),frameCnt/tBin]);
selAucMap = zeros([sum(~smallMask(:)),frameCnt/tBin]);
selErrorAucMap = zeros([sum(~smallMask(:)),frameCnt/tBin]);
aDPrime = zeros([sum(~smallMask(:)),frameCnt/tBin]);
cDPrime = zeros([sum(~smallMask(:)),frameCnt/tBin]);
selDPrime = zeros([sum(~smallMask(:)),frameCnt/tBin]);
selErrorDPrime = zeros([sum(~smallMask(:)),frameCnt/tBin]);

% compute dPrime
Cnt = 0;
for x = 1 : tBin : frameCnt
    
    Cnt = Cnt+1;
    temp = svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,x : x + tBin-1, 1:length(trials)),2))); %reconstruct frames
    temp = arrayShrink(arrayResize(temp,2),smallMask,'merge'); %collapse first 2 dimensions
    

    aAucMap(:,Cnt) = colAUC(temp(:,allInd)',aInd(allInd),'abs',false)' - 0.5;
    cAucMap(:,Cnt) = colAUC(temp(:,sucInd)',aInd(sucInd),'abs',false)' - 0.5;
    selAucMap(:,Cnt) = colAUC(temp(:,selInd)',aInd(selInd),'abs',false)' - 0.5;
    selErrorAucMap(:,Cnt) = colAUC(temp(:,selErrorInd)',aInd(selErrorInd),'abs',false)' - 0.5;

    aDPrime(:,Cnt) = (mean(temp(:,allInd & vInd),2) - mean(temp(:,allInd & aInd),2)) ./ std(temp(:,allInd), [],2);
    cDPrime(:,Cnt) = (mean(temp(:,sucInd & vInd),2) - mean(temp(:,sucInd & aInd),2)) ./ std(temp(:,sucInd), [],2);
    selDPrime(:,Cnt) = (mean(temp(:,selInd & vInd),2) - mean(temp(:,selInd & aInd),2)) ./ std(temp(:,selInd), [],2);
    selErrorDPrime(:,Cnt) = (mean(temp(:,selErrorInd & vInd),2) - mean(temp(:,selErrorInd & aInd),2)) ./ std(temp(:,selInd), [],2);
    
    disp(['Completed bin ' num2str(Cnt) '/' num2str(frameCnt/tBin)]);
end

save([cPath filesep Animal '-' Rec '-aDPrime.mat'],'aDPrime');
save([cPath filesep Animal '-' Rec '-cDPrime.mat'],'cDPrime');
save([cPath filesep Animal '-' Rec '-selDPrime.mat'],'selDPrime');
save([cPath filesep Animal '-' Rec '-selErrorDPrime.mat'],'selErrorDPrime');

save([cPath filesep Animal '-' Rec '-aAUC.mat'],'aAucMap');
save([cPath filesep Animal '-' Rec '-cAUC.mat'],'cAucMap');
save([cPath filesep Animal '-' Rec '-selAucMap.mat'],'selAucMap');
save([cPath filesep Animal '-' Rec '-selErrorAucMap.mat'],'selErrorAucMap');

%%
% allAuc(yy,:,:,:) = aucMap;
% save([cPath filesep Animal '_' Rec '_' figTitle{yy} '.mat'],'aucMap')
% smallMask = isnan(aucMap(:,:,1));
% aucMap = abs(arrayShrink(aucMap,smallMask,'merge')-0.5); %split pixels to get auc map
%
% numClust = 10;
% cMap = get(gca,'ColorOrder');
%
% Z = linkage(zscore(aucMap,[],2), 'ward','euclidean');
% c = cluster(Z,'maxclust',numClust);
% %         [~, T] = dendrogram(Z, numClust);
%
% figure('name',['Clust_' figTitle{yy}]);
% tI = 1:2:numClust*2;
% for x = 1:numClust
%     subplot(numClust,2,tI(x));
%     plot(mean(aucMap(c == x,:)));hold on; axis square
% end
% subplot(numClust,2,2:2:numClust*2);
% clustAuc = arrayShrink(c,smallMask,'split'); %split pixels to get auc map
% cImg = imshow(smooth2a(clustAuc,1));axis image;colormap jet;caxis([0 numClust])
% set(cImg,'AlphaData',~isnan(clustAuc)); %make NaNs transparent.
% set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
% saveas(gcf,[cPath filesep Animal '_' Rec '_Clust_' figTitle{yy} '.pdf']);
%
% figure('name',figTitle{yy})
% ind = frameCnt/4/tBin+1:frameCnt/size(allAuc,4)/tBin:frameCnt/tBin;
% for x = 1:length(ind)
%     subplot(3,size(allAuc,4)/4,x);
%     cImg = imshow(squeeze(allAuc(yy,:,:,ind(x))));axis image;caxis([0.15 0.75]); colormap jet
%     set(cImg,'AlphaData',~isnan(squeeze(allAuc(yy,:,:,ind(x))))); %make NaNs transparent.
%     title(['Time: ' num2str(time(ind(x)-1)) 's'])
% end
%
% set(gcf,'Position',[50 50 1200 800]);
% set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
% saveas(gcf,[cPath filesep Animal '_' Rec '_' figTitle{yy} '.pdf']);

%% AUC analysis
% cMap = get(gca,'ColorOrder');
% figTitle = {'AucVision' 'AucAudio' 'AucAudioVisual'};
% numClust = 10;
% allAuc = zeros([3,size(arrayResize(mask,4)),size(Vc,2)/tBin]);
%
% leftInd = SessionData.Assisted(trials) & ... %index for leftward CHOICES
%     ((SessionData.CorrectSide(trials) == 1 & SessionData.Rewarded(trials)) | ... %animal went left and got rewarded
%     (SessionData.CorrectSide(trials) == 2 & SessionData.Punished(trials))) ; %animal went left instead of right
%
% for yy = 1:3
%
%     ind = SessionData.StimType(trials) == yy; %find current modality
%     frameCnt = size(Vc,2);
%     time = ((1/sRate)*tBin:(1/sRate)*tBin:4)'-2;
%     classLoss = zeros(1,frameCnt/tBin);
%     cMap = get(gca,'ColorOrder');
%     smallMask = arrayResize(mask,4) == 1;
%     aucMap = zeros([sum(~smallMask(:)),frameCnt/tBin]);
%
%     Cnt = 0;
%     for x=1:tBin:frameCnt
%         Cnt = Cnt+1;
%         temp = svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,x:x+tBin-1,ind),2))); %reconstruct frames
%         temp = arrayShrink(arrayResize(temp,4),smallMask,'merge'); %collapse first 2 dimensions
%         aucMap(:,Cnt) = colAUC(temp',leftInd(ind),'abs',false)';
%         disp(['AUC map: Completed bin ' num2str(Cnt) '/' num2str(frameCnt/tBin)]);
%     end
%     aucMap = imrotate(arrayShrink(aucMap,smallMask,'split'),rotAngle,'crop');
%     aucMap(aucMap == 0) = NaN;
%     allAuc(yy,:,:,:) = aucMap;
%     save([cPath filesep Animal '_' Rec '_' figTitle{yy} '.mat'],'aucMap')
%
%     smallMask = isnan(aucMap(:,:,1));
%     aucMap = abs(arrayShrink(aucMap,smallMask,'merge')-0.5); %split pixels to get auc map
%
%     Z = linkage(zscore(aucMap,[],2), 'ward','euclidean');
%     c = cluster(Z,'maxclust',numClust);
% %         [~, T] = dendrogram(Z, numClust);
%
%     figure('name',['Clust_' figTitle{yy}]);
%     tI = 1:2:numClust*2;
%     for x = 1:numClust
%         subplot(numClust,2,tI(x));
%         plot(mean(aucMap(c == x,:)));hold on; axis square
%     end
%     subplot(numClust,2,2:2:numClust*2);
%     clustAuc = arrayShrink(c,smallMask,'split'); %split pixels to get auc map
%     cImg = imshow(smooth2a(clustAuc,1));axis image;colormap jet;caxis([0 numClust])
%     set(cImg,'AlphaData',~isnan(clustAuc)); %make NaNs transparent.
%     set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
%     saveas(gcf,[cPath filesep Animal '_' Rec '_Clust_' figTitle{yy} '.pdf']);
%
%     figure('name',figTitle{yy})
%     ind = frameCnt/4/tBin+1:frameCnt/size(allAuc,4)/tBin:frameCnt/tBin;
%     for x = 1:length(ind)
%         subplot(3,size(allAuc,4)/4,x);
%         cImg = imshow(squeeze(allAuc(yy,:,:,ind(x))));axis image;caxis([0.15 0.75]); colormap jet
%         set(cImg,'AlphaData',~isnan(squeeze(allAuc(yy,:,:,ind(x))))); %make NaNs transparent.
%         title(['Time: ' num2str(time(ind(x)-1)) 's'])
%     end
%
%     set(gcf,'Position',[50 50 1200 800]);
%     set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
%     saveas(gcf,[cPath filesep Animal '_' Rec '_' figTitle{yy} '.pdf']);
%
% end

%% build modality classifier
% h = figure('name',['SVM performance: ' Animal '-' Rec]); hold on;
% figTitle = {'Vision-Multisensory' 'Vision-Audio' 'Multisensory-Audio'};
% modBeta = zeros([3,size(mask),size(Vc,2)/tBin]);
%
% for xx = 1:3
%     if xx == 1
%         ind = SessionData.StimType(trials) ~= 2; %no audio trials
%         ind1 = SessionData.StimType(trials) == 1; %compare vision to multisensory trials
%     elseif xx == 2
%         ind = SessionData.StimType(trials) ~= 3; %no multisensory trials
%         ind1 = SessionData.StimType(trials) == 1; %compare vision to audio trials
%     elseif xx == 3
%         ind = SessionData.StimType(trials) ~= 1; %no vision trials
%         ind1 = SessionData.StimType(trials) == 3; %compare multisensory to audio trials
%     end
%
%     frameCnt = size(Vc,2);
%     time = ((1/sRate)*tBin:(1/sRate)*tBin:4)'-2;
%     betaMap = zeros([dimCnt,length(time)]);
%     classLoss = zeros(1,frameCnt/tBin);
%     cMap = get(gca,'ColorOrder');
%
%     Cnt = 0;
%     for x=1:tBin:frameCnt
%         Cnt = Cnt+1;
%
%         SVMModel = fitcsvm(squeeze(mean(Vc(1:dimCnt,x:x+tBin-1,ind),2))',ind1(ind),'KernelFunction','linear','Standardize',true,'ClassNames',[false true]);
%         CVSVMModel = crossval(SVMModel,'kfold',crossValCnt); %cross-validate model
%         classLoss(Cnt) = kfoldLoss(CVSVMModel); %compute prediction error
%
%         temp = zeros([crossValCnt,dimCnt]);
%         for iModels = 1:crossValCnt
%             temp(iModels,:) = CVSVMModel.Trained{iModels}.Beta; %get beta weights from individual training models
%         end
%         betaMap(1:dimCnt,Cnt) = mean(temp); %average over beta weights from all trained models
%         disp(['SVM map: Completed bin ' num2str(Cnt) '/' num2str(frameCnt/tBin)]);
%
%     end
%     save([cPath filesep Animal '_' Rec '_' figTitle{xx} '.mat'],'SVMModel','CVSVMModel','classLoss','betaMap')
%
%     plot(h.Children,time,smooth(1-classLoss));
%     modBeta(xx,:,:,:) = imrotate(svdFrameReconstruct(U(:,:,1:dimCnt),betaMap),rotAngle,'crop');
%
%     figure('name',figTitle{xx})
%     ind = [frameCnt/2/tBin+1:frameCnt/size(modBeta,4)/tBin:frameCnt/tBin frameCnt/tBin];
%     for x = 1:size(modBeta,4)/2
%         subplot(2,size(modBeta,4)/4,x);
%         imagesc(smooth2a(squeeze(abs(modBeta(xx,:,:,x))),2));axis image;caxis([0 0.02]); colormap jet
%         title(['Time: ' num2str(time(ind(x))) 's'])
%     end
%     set(gcf,'Position',[50 50 1200 800]);
%     set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
%     saveas(gcf,[cPath filesep Animal '_' figTitle{xx} '.pdf']);
%
% end
% modBeta(modBeta == 0) = NaN;
%
% figure(h);
% legend(figTitle);axis square; ylim([0.4 1]);
% title(['SVM performance: ' Animal '-' Rec])
% xlabel('Time (s)'); ylabel('Model Performance (%)')
% title(['SVM performance: ' Animal '-' Rec]);
% set(gcf,'Position',[50 50 1200 800]);
% set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
% saveas(h,[cPath filesep 'Modalityperformance-' Animal '-' Rec '-' figTitle{xx} '.pdf']);
% save([cPath filesep 'ModBetaMatrix-' Animal '-' Rec '.mat'],'modBeta');
%
% %% Decision classifier
% h = figure('name',['SVM performance: ' Animal '-' Rec]); hold on;
% cMap = get(gca,'ColorOrder');
% figTitle = {'ChoiceVision' 'ChoiceAudio' 'ChoiceAudioVisual'};
% choiceBeta = zeros([3,size(mask),size(Vc,2)/tBin]);
%
% f = @(p,x) p(4) + p(1) ./ (1 + exp(-(x-p(2))/p(3)));
% xx = linspace(-2,2)';
%
% leftInd = SessionData.Assisted(trials) & ... %index for leftward CHOICES
%     ((SessionData.CorrectSide(trials) == 1 & SessionData.Rewarded(trials)) | ... %animal went left and got rewarded
%     (SessionData.CorrectSide(trials) == 2 & SessionData.Punished(trials))) ; %animal went left instead of right
%
% for yy = 1:3
%
%     ind = SessionData.StimType(trials) == yy; %find current modality
%     frameCnt = size(Vc,2);
%     time = ((1/sRate)*tBin:(1/sRate)*tBin:4)'-2;
%     betaMap = zeros([dimCnt,length(time)]);
%     classLoss = zeros(1,frameCnt/tBin);
%     cMap = get(gca,'ColorOrder');
%
%     Cnt = 0;
%     for x=1:tBin:frameCnt
%         Cnt = Cnt+1;
%
%         SVMModel = fitcsvm(squeeze(mean(Vc(1:dimCnt,x:x+tBin-1,ind),2))',leftInd(ind),'KernelFunction','linear','Standardize',true,'ClassNames',[false true]);
%         CVSVMModel = crossval(SVMModel,'kfold',crossValCnt); %cross-validate model
%         classLoss(Cnt) = kfoldLoss(CVSVMModel); %compute prediction error
%
%         temp = zeros([crossValCnt,dimCnt]);
%         for iModels = 1:crossValCnt
%             temp(iModels,:) = CVSVMModel.Trained{iModels}.Beta; %get beta weights from individual training models
%         end
%         betaMap(1:dimCnt,Cnt) = mean(temp); %average over beta weights from all trained models
%         disp(['SVM map: Completed bin ' num2str(Cnt) '/' num2str(frameCnt/tBin)]);
%
%     end
%     save([cPath filesep Animal '_' Rec '_' figTitle{yy} '.mat'],'SVMModel','CVSVMModel','classLoss','betaMap')
%
%     %plot data
%     plot(h.Children,time,smooth(1-classLoss));
% %     plot(h.Children,time,smooth(1-classLoss),'o','Markersize',5,'MarkerFaceColor','w','linewidth',2,'color',cMap(yy,:));
% %     Est = [classLoss(1) time(round(length(time)*0.6)) (max(time)-min(time))/10 0.5];
% %     line(h.Children,xx,1-predict(fitnlm(time,classLoss,f,Est),xx),'linewidth',2,'color',cMap(yy,:)); % plot fit
%
%     choiceBeta(yy,:,:,:) = imrotate(svdFrameReconstruct(U(:,:,1:dimCnt),betaMap),rotAngle,'crop');
%
%     figure('name',figTitle{yy})
%     ind = [frameCnt/2/tBin+1:frameCnt/size(choiceBeta,4)/tBin:frameCnt/tBin frameCnt/tBin];
%     for x = 1:size(choiceBeta,4)/2
%         subplot(2,size(choiceBeta,4)/4,x);
%         imagesc(smooth2a(squeeze(abs(choiceBeta(yy,:,:,x))),2));axis image;caxis([0 0.03]); colormap jet
%         title(['Time: ' num2str(time(ind(x))) 's'])
%     end
%
% %     ind = [frameCnt/2/tBin+1:frameCnt/16/tBin:frameCnt/tBin frameCnt/tBin];
% %     for x = 1:8
% %         subplot(2,4,x);
% %         imagesc(smooth2a(mean(abs(temp(:,:,ind(x):ind(x+1))),3),2));axis image;caxis([0 0.005]); colormap jet
% %         title(['Time: ' num2str(time(ind(x))) 's'])
% %     end
%
%     set(gcf,'Position',[50 50 1200 800]);
%     set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
%     saveas(gcf,[cPath filesep Animal '_' Rec '_' figTitle{yy} '.pdf']);
% end
%
% choiceBeta(choiceBeta == 0) = NaN;
%
% figure(h);
% legend(figTitle);axis square; ylim([0.4 1]);
% title(['SVM performance: ' Animal '-' Rec])
% xlabel('Time (s)'); ylabel('Model Performance (%)')
% title(['SVM performance: ' Animal '-' Rec]);
% set(gcf,'Position',[50 50 1200 800]);
% set(gcf,'PaperOrientation','landscape','PaperPositionMode','auto');
% saveas(h,[cPath filesep 'ChoicePerformance_' Animal '_' Rec '.pdf']);
% save([cPath filesep 'ChoiceBetaMatrix-' Animal '-' Rec '.mat'],'choiceBeta');


%% betaMapSingle(:,:,Cnt,iDist) = Widefield_DimMerge_v1(SVMModel.Beta,mask,'split'); %split pixels into original frame size
%
%
% for iDist = 1:length(distStim)
%     Cnt = 0;
%     for x=1:tBin:frameCnt
%         Cnt = Cnt+1;
%
%         detectInd = SessionData.StimType(trials) == 3 & distFractions(trials) == distStim(iDist); %multisensory trials at current distFrequency
%         SVMModel = fitcsvm(squeeze(mean(binData(:,x:x+tBin-1,detectInd),2))',leftInd(detectInd),'KernelFunction','linear','Standardize',true,'ClassNames',{'false','true'});
% %         betaMap(:,:,Cnt,iDist) = Widefield_DimMerge_v1(SVMModel.Beta,mask,'split'); %split pixels into original frame size
%
%        CVSVMModel = crossval(SVMModel); %cross-validate mode;
%        classLoss(Cnt,iDist) = kfoldLoss(CVSVMModel); %compute prediction error
%
%        disp(['SVM map: Completed bin ' num2str(Cnt) '/' num2str(frameCnt/tBin) '; DistFractions = ' num2str(distStim(iDist))]);
%         %     SVMModel = fitclinear(squeeze(mean(binData(:,x:x+tBin-1,:),2)),iLeft,'ClassNames',{'false','true'},'ObserVationsIn','columns');
%     end
%     Est = [classLoss(1,iDist) time(round(length(time)*0.6)) (max(time)-min(time))/10 0.5];
%     linearFit{iDist} = fitnlm(time,classLoss(:,iDist),f,Est);
%
%     %plot data
%     subplot(1,length(distStim),iDist);
%     plot(time,1-classLoss(:,iDist),'o','Markersize',10,'MarkerFaceColor','w','linewidth',2,'color',cMap(iDist,:));
%     line(xx,1-predict(linearFit{iDist},xx),'linewidth',2,'color',cMap(iDist,:)); % plot fit
%     ylim([0.4 1]);xlim([-1 2]);title(['Dist. fraction: ' num2str(distStim(iDist))]);
%     xlabel('time (s)');ylabel('Prediction error');axis square
%     temp = sort(1-predict(linearFit{iDist},xx)); %sort clasifier performance
%     modelPerf(iDist) = median(temp(end-4:end));
%
%     estAmps(iDist) = linearFit{iDist}.Coefficients.Estimate(2); %M50 - inflection point
%     line([estAmps(iDist),estAmps(iDist)],[-1 2],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
% end
%
% figure(bhvFig)
% subplot(1,2,1);
% plot(distStim,modelPerf,'-o','Color',[.5 .5 .5],'Markersize',9,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor','w','linewidth',2)
% subplot(1,2,2);
% plot(distStim,estAmps,'-o','Color',[.5 .5 .5],'Markersize',9,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor','w','linewidth',2)
% title([Animal ' - Clasifier inflection point time']); xlabel('Dist. fraction','fontsize',15); ylabel('time(s)','fontsize',15);
% axis square; set(gca,'FontSize',12); xlim([-.05 .85]);ylim([0 2]);
%
%
% %% build choice classifier - based on all trials were animal chose left vs right
% leftInd = SessionData.Assisted(trials) & ...
%     ((SessionData.CorrectSide(trials) == 1 & SessionData.Rewarded(trials)) | ... %animal went left and got rewarded
%     (SessionData.CorrectSide(trials) == 2 & SessionData.Punished(trials))) ; %animal went left instead of right
%
% tBin = 1;
% time = ((4/trialDur)*tBin:(4/trialDur)*tBin:4)'-2;
% f = @(p,x) p(4) + p(1) ./ (1 + exp(-(x-p(2))/p(3)));
% xx = linspace(-1,2)';
%
% trialDur = size(binData,2);
% % betaMap = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
% classLoss = zeros(trialDur/tBin,length(distStim));
% figure;axis square;hold on;
% cMap = get(gca,'ColorOrder');
%
% for iDist = 1:length(distStim)
%     Cnt = 0;
%     for x=1:tBin:trialDur
%         Cnt = Cnt+1;
%
%         detectInd = SessionData.StimType(trials) == 3 & distFractions(trials) == distStim(iDist); %multisensory trials at current distFrequency
%         SVMModel = fitcsvm(squeeze(mean(binData(:,x:x+tBin-1,detectInd),2))',leftInd(detectInd),'KernelFunction','linear','Standardize',true,'ClassNames',{'false','true'});
% %         betaMap(:,:,Cnt,iDist) = Widefield_DimMerge_v1(SVMModel.Beta,mask,'split'); %split pixels into original frame size
%
%         CVSVMModel = crossval(SVMModel); %cross-validate mode;
%         classLoss(Cnt,iDist) = kfoldLoss(CVSVMModel); %compute prediction error
%
%         disp(['SVM map: Completed bin ' num2str(Cnt) '/' num2str(trialDur/tBin) '; DistFractions = ' num2str(distStim(iDist))]);
%         %     SVMModel = fitclinear(squeeze(mean(binData(:,x:x+tBin-1,:),2)),iLeft,'ClassNames',{'false','true'},'ObserVationsIn','columns');
%     end
%     Est = [classLoss(1,iDist) time(round(length(time)*0.75)) (max(time)-min(time))/10 0.5];
%     linearFit{iDist} = fitnlm(time,classLoss(:,iDist),f,Est);
%
%     %plot data
%     subplot(1,length(distStim),iDist);
%     plot(time,1-classLoss(:,iDist),'o','Markersize',10,'MarkerFaceColor','w','linewidth',2,'color',cMap(iDist,:));
%     line(xx,1-predict(linearFit{iDist},xx),'linewidth',2,'color',cMap(iDist,:)); % plot fit
%     ylim([0.4 1]);xlim([-1 2]);title(['Dist. fraction: ' num2str(distStim(iDist))]);
%     xlabel('time (s)');ylabel('Prediction error');axis square
%
%     p = linearFit{iDist}.Coefficients.Estimate';
%     estAmps(iDist) = p(2);                  %M50 - inflection point
%     line([estAmps(iDist),estAmps(iDist)],[-1 2],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
% end
%
%
% %% build choice classifier - based on all trials were animal correctly chose left vs right
% leftInd = SessionData.Assisted(trials) & ...
%     ((SessionData.CorrectSide(trials) == 1 & SessionData.Rewarded(trials)) | ... %animal went left and got rewarded
%     (SessionData.CorrectSide(trials) == 2 & SessionData.Punished(trials))) ; %animal went left instead of right
%
% tBin = 1;
% time = ((4/trialDur)*tBin:(4/trialDur)*tBin:4)'-2;
% f = @(p,x) p(4) + p(1) ./ (1 + exp(-(x-p(2))/p(3)));
% xx = linspace(-1,2)';
%
% trialDur = size(binData,2);
% % betaMap = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
% classLossTrain = zeros(trialDur/tBin,length(distStim));
% classLossPredict = zeros(trialDur/tBin,length(distStim));
% figure;axis square;
% cMap = get(gca,'ColorOrder');
%
% for iDist = 1:length(distStim)
%     Cnt = 0;
%     for x=1:tBin:trialDur
%         Cnt = Cnt+1;
%
%         correctInd = SessionData.StimType(trials) == 3 & distFractions(trials) == distStim(iDist) & SessionData.Rewarded(trials); %multisensory trials at current distFrequency. Train on correct choice only
%         trialCountTrain(Cnt,iDist) = sum(correctInd);
%         SVMModel = fitcsvm(squeeze(mean(binData(:,x:x+tBin-1,correctInd),2))',leftInd(correctInd),'KernelFunction','linear','Standardize',true,'ClassNames',{'false','true'});
%         %         betaMap(:,:,Cnt,iDist) = Widefield_DimMerge_v1(SVMModel.Beta,mask,'split'); %split pixels into original frame size
%
%         CVSVMModel = crossval(SVMModel); %cross-validate mode;
%         classLossTrain(Cnt,iDist) = kfoldLoss(CVSVMModel); %compute prediction error
%
%         errorInd = SessionData.StimType(trials) == 3 & distFractions(trials) == distStim(iDist) & ~SessionData.Rewarded(trials); %multisensory trials at current distFrequency. Test on error choice only
%         trialCountPredict(Cnt,iDist) = sum(errorInd);
%         classLossPredict(Cnt,iDist) = mean((ismember(predict(SVMModel,squeeze(mean(binData(:,x:x+tBin-1,errorInd),2))'),'true') - leftInd(errorInd)').^2); %try to predict outcome in error trials
%
%         disp(['SVM map: Completed bin ' num2str(Cnt) '/' num2str(trialDur/tBin) '; DistFractions = ' num2str(distStim(iDist))]);
%         %     SVMModel = fitclinear(squeeze(mean(binData(:,x:x+tBin-1,:),2)),iLeft,'ClassNames',{'false','true'},'ObserVationsIn','columns');
%     end
%     Est = [classLossTrain(1,iDist) time(round(length(time)*0.7)) (max(time)-min(time))/50 0.5];
%     weights = ones(1,length(time));
% %     weights(time>1 & time<=1.5) = 50;
%     linearFitTrain{iDist} = fitnlm(time,classLossTrain(:,iDist),f,Est,'Weights',weights);
%
%     %plot trained data
%     subplot(2,length(distStim),iDist);
%     plot(time,1-classLossTrain(:,iDist),'o','Markersize',10,'MarkerFaceColor','w','linewidth',2,'color',cMap(iDist,:));
%     line(xx,1-predict(linearFitTrain{iDist},xx),'linewidth',2,'color',cMap(iDist,:)); % plot fit
%     ylim([0.4 1]);xlim([-1 2]);title(['Dist. fraction: ' num2str(distStim(iDist))]);
%     xlabel('time (s)');ylabel('Prediction error');axis square
%
%     p = linearFitTrain{iDist}.Coefficients.Estimate';
%     estAmps{iDist,1}(1) = p(2)-((p(3)*10)/2);    %M0 - onset of sigmoid rise
%     estAmps{iDist,1}(2) = p(2)-((p(3)*2.2)/2);   %M25 - 25% of sigmoid rise
%     estAmps{iDist,1}(3) = p(2);                  %M50 - inflection point
%     estAmps{iDist,1}(4) = ((p(3)*2.2)/2)+p(2);   %M75 75% of sigmoid rise
%     estAmps{iDist,1}(5) = ((p(3)*10)/2)+p(2);    %M100 - onset of sigmoid plateu
%     estAmps{iDist,1}(estAmps{iDist,1} <= 0) = 0;  % amps can't be below 0
%     line([estAmps{iDist,1}(3),estAmps{iDist,1}(3)],[-1 2],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
%
%     %plot predicted data
%     Est = [classLossPredict(1,iDist) time(round(length(time)*0.75)) (max(time)-min(time))/10 0.5];
%     linearFitPredict{iDist} = fitnlm(time,classLossPredict(:,iDist),f,Est);
%
%     subplot(2,length(distStim),iDist+length(distStim));
%     plot(time,1-classLossPredict(:,iDist),'o','Markersize',10,'MarkerFaceColor','w','linewidth',2,'color',cMap(iDist,:));
%     line(xx,1-predict(linearFitPredict{iDist},xx),'linewidth',2,'color',cMap(iDist,:)); % plot fit
%     ylim([0.4 1]);xlim([-1 2]);title(['Dist. fraction: ' num2str(distStim(iDist))]);
%     xlabel('time (s)');ylabel('Prediction error');axis square
%
%     p = linearFitPredict{iDist}.Coefficients.Estimate';
%     estAmps{iDist,2}(1) = p(2)-((p(3)*10)/2);    %M0 - onset of sigmoid rise
%     estAmps{iDist,2}(2) = p(2)-((p(3)*2.2)/2);   %M25 - 25% of sigmoid rise
%     estAmps{iDist,2}(3) = p(2);                  %M50 - inflection point
%     estAmps{iDist,2}(4) = ((p(3)*2.2)/2)+p(2);   %M75 75% of sigmoid rise
%     estAmps{iDist,2}(5) = ((p(3)*10)/2)+p(2);    %M100 - onset of sigmoid plateu
%     estAmps{iDist,2}(estAmps{iDist,1} <= 0) = 0;  % amps can't be below 0
%     line([estAmps{iDist,2}(3),estAmps{iDist,2}(3)],[-1 2],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
% end
%
%
% %% build choice classifier - perform random crossvalidation 100 times and average beta maps.
% leftInd = SessionData.Assisted(trials) & ...
%     ((SessionData.CorrectSide(trials) == 1 & SessionData.Rewarded(trials)) | ... %animal went left and got rewarded
%     (SessionData.CorrectSide(trials) == 2 & SessionData.Punished(trials))) ; %animal went left instead of right
% % leftInd = SessionData.Assisted(trials) & SessionData.CorrectSide(trials) == 1; %leftward trials
% tBin = 1;
% crossValCnt = 10; %amount of repetitions for beta maps
% holdOut = 0.5; %percent of data that should be used for crossvalidation
%
% trialDur = size(binData,2);
% betaMap = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
% betaStd = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
% classLoss = zeros(trialDur/tBin,length(distStim));
% figure;axis square;hold on;
% cMap = get(gca,'ColorOrder');
%
% % for iDist = 1:length(distStim)
% for iDist = length(distStim)
%     Cnt = 0;
%     for x=1:tBin:trialDur
%         Cnt = Cnt+1;
%
%         temp = zeros(size(binData,1),crossValCnt);
%         for iModel = 1:crossValCnt
%
%             trainSelect = rand(1,sum(trials)) > holdOut; %randomly select subset of all trials for training
%             detectInd = SessionData.StimType(trials) == 3 & distFractions(trials) >= distStim(iDist) & trainSelect; %multisensory trials at current distFrequency
%             SVMModel = fitcsvm(squeeze(mean(binData(:,x:x+tBin-1,detectInd),2))',leftInd(detectInd),'KernelFunction','linear','Standardize',true);
%             temp(:,iModel) = SVMModel.Beta;
%             testInd = SessionData.StimType(trials) == 3 & distFractions(trials) == distStim(iDist) & ~trainSelect; %use prediction subset to test current model
%             modelLoss(iModel) = mean((leftInd(testInd)'-predict(SVMModel,squeeze(mean(binData(:,x:x+tBin-1,testInd),2))')).^2);
%
%         end
%
%         betaMapSingle(:,:,Cnt,iDist) = Widefield_DimMerge_v1(SVMModel.Beta,mask,'split'); %split pixels into original frame size
%         betaMap(:,:,Cnt,iDist) = Widefield_DimMerge_v1(mean(temp,2),mask,'split'); %split pixels into original frame size
%         betaStd(:,:,Cnt,iDist) = Widefield_DimMerge_v1(std(temp,0,2),mask,'split'); %split pixels into original frame size
%         classLoss(Cnt,iDist) = mean(modelLoss); %compute prediction error
%
%         disp(['SVM map: Completed bin ' num2str(Cnt) '/' num2str(trialDur/tBin) '; DistFractions = ' num2str(distStim(iDist))]);
%         %     SVMModel = fitclinear(squeeze(mean(binData(:,x:x+tBin-1,:),2)),iLeft,'ClassNames',{'false','true'},'ObserVationsIn','columns');
%     end
% end
%
% %% plot beta maps
% % combinedBeta = abs(betaMap./betaStd); %normalize betaMap with standard deviaton
%
% for iDist = 1:5
%     figure
%     subplot(1,3,1)
%     imagesc(smooth2a(mean(abs(betaMap(:,:,:,iDist)),3),2));axis square;colormap jet
%     title(['Distractor: ' num2str(distStim(iDist)) '; All frames'])
%     subplot(1,3,2)
%     imagesc(smooth2a(mean(abs(betaMap(:,:,80:120,iDist)),3),2));axis square;colormap jet
%     title(['Distractor: ' num2str(distStim(iDist)) '; Frame 40-50'])
%     subplot(1,3,3)
%     imagesc(smooth2a(mean(abs(betaMap(:,:,120:160,iDist)),3),2));axis square;colormap jet
%     title(['Distractor: ' num2str(distStim(iDist)) '; Frame 60-80'])
% end
%
%
% %% compute auc maps
% %% build choice classifier - perform random crossvalidation 100 times and average beta maps.
% leftInd = SessionData.Assisted(trials) & ...
%     ((SessionData.CorrectSide(trials) == 1 & SessionData.Rewarded(trials)) | ... %animal went left and got rewarded
%     (SessionData.CorrectSide(trials) == 2 & SessionData.Punished(trials))) ; %animal went left instead of right
% % leftInd = SessionData.Assisted(trials) & SessionData.CorrectSide(trials) == 1; %leftward trials
% tBin = 1;
% crossValCnt = 100; %amount of repetitions for beta maps
% holdOut = 0.5; %percent of data that should be used for crossvalidation
%
% trialDur = size(binData,2);
% aucMap = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
% aucMapStd = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
%
% % for iDist = 1:length(distStim)
% for iDist = [1 length(distStim)]
%     Cnt = 0;
%     for x=1:tBin:trialDur
%         Cnt = Cnt+1;
%
%         temp = zeros(size(binData,1),crossValCnt);
%         for iModel = 1:crossValCnt
%             trainSelect = rand(1,sum(trials)) > holdOut; %randomly select subset of all trials for training
%             detectInd = SessionData.StimType(trials) == 3 & distFractions(trials) == distStim(iDist) & trainSelect; %multisensory trials at current distFrequency
%             temp(:,iModel) = colAUC((squeeze(mean(binData(:,x:x+tBin-1,detectInd),2))./4095)',leftInd(detectInd),'abs',false);
%         if rem(iModel,10) == 0; disp(iModel);end
%         end
%
%         aucMap(:,:,Cnt,iDist) = Widefield_DimMerge_v1(mean(temp,2),mask,'split'); %split pixels into original frame size
%         aucMapStd(:,:,Cnt,iDist) = Widefield_DimMerge_v1(std(temp,0,2),mask,'split'); %split pixels into original frame size
%
%         disp(['AUC map: Completed bin ' num2str(Cnt) '/' num2str(trialDur/tBin) '; DistFractions = ' num2str(distStim(iDist))]);
%         %     SVMModel = fitclinear(squeeze(mean(binData(:,x:x+tBin-1,:),2)),iLeft,'ClassNames',{'false','true'},'ObserVationsIn','columns');
%     end
% end
% combinedAuc = abs(aucMap./aucMapStd); %normalize aucmap with standard deviaton