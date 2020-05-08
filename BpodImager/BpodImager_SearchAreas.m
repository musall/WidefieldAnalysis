function aInd = BpodImager_SearchAreas(Animal,Paradigm,Rec)
% Code to analyze stacks that are created by WidefieldMapperAnalyze_v3 and later.

%% basic variables
binSize = 4;
path = 'H:\BpodImager\Animals\'; %Widefield data path
smth = 1;           %2-D smooth for plots
rThresh = 98;       %treshold for ROI
cRange = [];        %color range for selection image. Leave empty to scale.1
stimOn = 10;        %frame at which stimulus was presented
predictAnimal = false;

%% load data
load([path Animal '\' Paradigm '\' Rec '\binData_' num2str(binSize) '.mat'])

%% load behavior data and get indices for different trial types
load([path Animal '\' Paradigm '\' Rec '\perfTrials.mat'])
cFile = ls([path Animal '\' Paradigm '\' Rec '\' Animal '_' Paradigm '*.mat']);
load([path Animal '\' Paradigm '\' Rec '\' cFile]);

if length(perfTrials) ~= length(SessionData.Rewarded)
    warning(['Inconsistent trialcount in perfTrials(' num2str(length(perfTrials)) ') and behavioral data(' num2str(length(SessionData.Rewarded)) '). Cut down for equal lengths.'])
    a = min([length(perfTrials) length(SessionData.Rewarded)]);
    perfTrials(a+1:end) = []; 
    SessionData.Assisted(a+1:end) = []; SessionData.Rewarded(a+1:end) = []; SessionData.CorrectSide(a+1:end) = [];
end

%% compute delta F/F movies
distFractions = SessionData.DistStim ./ SessionData.TargStim;
distStim = unique(distFractions);
leftInd = SessionData.Assisted(perfTrials) & SessionData.CorrectSide(perfTrials) == 1; %leftward trials

for iDist = 1:length(distStim)
    detectInd = SessionData.StimType(perfTrials) == 3 & distFractions(perfTrials) == distStim(iDist) & ~SessionData.DidNotChoose(perfTrials); %multisensory trials at current distFrequency
    
    % get all leftward trials, normalize against baseline and save as movie
    trace = Widefield_DimMerge_v1(mean(binData(:,:,leftInd(detectInd)),3),mask,'split');
    baseline = mean(trace(:,:,1:20),3);
    for iFrames = 1:size(trace,3)
        trace(:,:,iFrames) = (trace(:,:,iFrames) - baseline) ./ baseline;
    end
    % Widefield_SaveToAvi(trace,'LeftDetection',5,'jet',[0 double(max(trace(:)))/10],[path Animal '\' Paradigm '\'  Rec],1);
    Widefield_SaveToAvi(trace,['LeftDetection_' num2str(distStim(iDist))],5,'jet',[-0.02 0.02],[path Animal '\' Paradigm '\'  Rec],1);

    % get all rightward trials, normalize against baseline and save as movie
    trace = Widefield_DimMerge_v1(mean(binData(:,:,~leftInd(detectInd)),3),mask,'split');
    baseline = mean(trace(:,:,1:20),3);
    for iFrames = 1:size(trace,3)
        trace(:,:,iFrames) = (trace(:,:,iFrames) - baseline) ./ baseline;
    end
    Widefield_SaveToAvi(trace,['RightDetection_' num2str(distStim(iDist))],5,'jet',[-0.02 0.02],[path Animal '\' Paradigm '\'  Rec],1);
    clear trace baseline
end

%% compute animal performance and plot
for iDist = 1:length(distStim)
    ind = distFractions == distStim(iDist) & SessionData.StimType == 3; %only mulitsensory trials for initial performance curves
    [bhvPerf(iDist,1),bhvError(:,iDist,1),tCount(iDist,1)] = computeBehavior(SessionData,500,ind); %compute performance and error for all trials
end

bhvError(1,:,:) = bhvPerf - squeeze(bhvError(1,:,:))';
bhvError(2,:,:) = squeeze(bhvError(2,:,:))' - bhvPerf;
    
bhvFig = figure;
subplot(1,2,1);hold on
line([-0.2 1],[0.5 0.5],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
errorbar(distStim,bhvPerf,bhvError(1,:),bhvError(2,:), ...
    '-ok','Markersize',9,'MarkerEdgeColor','k','MarkerFaceColor','w','linewidth',2); 

title([Animal ' - Discrimination Performance']); xlabel(['Trialcounts: ' num2str(tCount')],'fontsize',15); ylabel('performance (%)','fontsize',15);
axis square; ylim([0.35 1]); xlim([-0.1 1.1]); set(gca,'FontSize',12); 

%% build side classifier
leftInd = SessionData.Assisted(perfTrials) & SessionData.CorrectSide(perfTrials) == 1; %leftward trials
trialDur = size(binData,2);
tBin = 1;
time = ((4/trialDur)*tBin:(4/trialDur)*tBin:4)'-2;
f = @(p,x) p(4) + p(1) ./ (1 + exp(-(x-p(2))/p(3)));
xx = linspace(-1,2)';

% betaMap = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
classLoss = zeros(trialDur/tBin,length(distStim));
figure;axis square;hold on;
cMap = get(gca,'ColorOrder');

for iDist = 1:length(distStim)
    Cnt = 0;
    for x=1:tBin:trialDur
        Cnt = Cnt+1;

        detectInd = SessionData.StimType(perfTrials) == 3 & distFractions(perfTrials) == distStim(iDist); %multisensory trials at current distFrequency
        SVMModel = fitcsvm(squeeze(mean(binData(:,x:x+tBin-1,detectInd),2))',leftInd(detectInd),'KernelFunction','linear','Standardize',true,'ClassNames',{'false','true'});
%         betaMap(:,:,Cnt,iDist) = Widefield_DimMerge_v1(SVMModel.Beta,mask,'split'); %split pixels into original frame size

       CVSVMModel = crossval(SVMModel); %cross-validate mode;
       classLoss(Cnt,iDist) = kfoldLoss(CVSVMModel); %compute prediction error

       disp(['SVM map: Completed bin ' num2str(Cnt) '/' num2str(trialDur/tBin) '; DistFractions = ' num2str(distStim(iDist))]);
        %     SVMModel = fitclinear(squeeze(mean(binData(:,x:x+tBin-1,:),2)),iLeft,'ClassNames',{'false','true'},'ObserVationsIn','columns');
    end
    Est = [classLoss(1,iDist) time(round(length(time)*0.6)) (max(time)-min(time))/10 0.5];
    linearFit{iDist} = fitnlm(time,classLoss(:,iDist),f,Est);
    
    %plot data
    subplot(1,length(distStim),iDist);
    plot(time,1-classLoss(:,iDist),'o','Markersize',10,'MarkerFaceColor','w','linewidth',2,'color',cMap(iDist,:));
    line(xx,1-predict(linearFit{iDist},xx),'linewidth',2,'color',cMap(iDist,:)); % plot fit
    ylim([0.4 1]);xlim([-1 2]);title(['Dist. fraction: ' num2str(distStim(iDist))]);
    xlabel('time (s)');ylabel('Prediction error');axis square
    temp = sort(1-predict(linearFit{iDist},xx)); %sort clasifier performance
    modelPerf(iDist) = median(temp(end-4:end));    
    
    estAmps(iDist) = linearFit{iDist}.Coefficients.Estimate(2); %M50 - inflection point
    line([estAmps(iDist),estAmps(iDist)],[-1 2],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
end

figure(bhvFig)
subplot(1,2,1);
plot(distStim,modelPerf,'-o','Color',[.5 .5 .5],'Markersize',9,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor','w','linewidth',2)
subplot(1,2,2);
plot(distStim,estAmps,'-o','Color',[.5 .5 .5],'Markersize',9,'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor','w','linewidth',2)
title([Animal ' - Clasifier inflection point time']); xlabel('Dist. fraction','fontsize',15); ylabel('time(s)','fontsize',15);
axis square; set(gca,'FontSize',12); xlim([-.05 .85]);ylim([0 2]);

       
%% build choice classifier - based on all trials were animal chose left vs right
leftInd = SessionData.Assisted(perfTrials) & ...
    ((SessionData.CorrectSide(perfTrials) == 1 & SessionData.Rewarded(perfTrials)) | ... %animal went left and got rewarded
    (SessionData.CorrectSide(perfTrials) == 2 & SessionData.Punished(perfTrials))) ; %animal went left instead of right

tBin = 1;
time = ((4/trialDur)*tBin:(4/trialDur)*tBin:4)'-2;
f = @(p,x) p(4) + p(1) ./ (1 + exp(-(x-p(2))/p(3)));
xx = linspace(-1,2)';

trialDur = size(binData,2);
% betaMap = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
classLoss = zeros(trialDur/tBin,length(distStim));
figure;axis square;hold on;
cMap = get(gca,'ColorOrder');

for iDist = 1:length(distStim)
    Cnt = 0;
    for x=1:tBin:trialDur
        Cnt = Cnt+1;
        
        detectInd = SessionData.StimType(perfTrials) == 3 & distFractions(perfTrials) == distStim(iDist); %multisensory trials at current distFrequency
        SVMModel = fitcsvm(squeeze(mean(binData(:,x:x+tBin-1,detectInd),2))',leftInd(detectInd),'KernelFunction','linear','Standardize',true,'ClassNames',{'false','true'});
%         betaMap(:,:,Cnt,iDist) = Widefield_DimMerge_v1(SVMModel.Beta,mask,'split'); %split pixels into original frame size
        
        CVSVMModel = crossval(SVMModel); %cross-validate mode;
        classLoss(Cnt,iDist) = kfoldLoss(CVSVMModel); %compute prediction error
        
        disp(['SVM map: Completed bin ' num2str(Cnt) '/' num2str(trialDur/tBin) '; DistFractions = ' num2str(distStim(iDist))]);
        %     SVMModel = fitclinear(squeeze(mean(binData(:,x:x+tBin-1,:),2)),iLeft,'ClassNames',{'false','true'},'ObserVationsIn','columns');
    end
    Est = [classLoss(1,iDist) time(round(length(time)*0.75)) (max(time)-min(time))/10 0.5];
    linearFit{iDist} = fitnlm(time,classLoss(:,iDist),f,Est);
    
    %plot data
    subplot(1,length(distStim),iDist);
    plot(time,1-classLoss(:,iDist),'o','Markersize',10,'MarkerFaceColor','w','linewidth',2,'color',cMap(iDist,:));
    line(xx,1-predict(linearFit{iDist},xx),'linewidth',2,'color',cMap(iDist,:)); % plot fit
    ylim([0.4 1]);xlim([-1 2]);title(['Dist. fraction: ' num2str(distStim(iDist))]);
    xlabel('time (s)');ylabel('Prediction error');axis square
    
    p = linearFit{iDist}.Coefficients.Estimate';
    estAmps(iDist) = p(2);                  %M50 - inflection point
    line([estAmps(iDist),estAmps(iDist)],[-1 2],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
end

       
%% build choice classifier - based on all trials were animal correctly chose left vs right
leftInd = SessionData.Assisted(perfTrials) & ...
    ((SessionData.CorrectSide(perfTrials) == 1 & SessionData.Rewarded(perfTrials)) | ... %animal went left and got rewarded
    (SessionData.CorrectSide(perfTrials) == 2 & SessionData.Punished(perfTrials))) ; %animal went left instead of right

tBin = 1;
time = ((4/trialDur)*tBin:(4/trialDur)*tBin:4)'-2;
f = @(p,x) p(4) + p(1) ./ (1 + exp(-(x-p(2))/p(3)));
xx = linspace(-1,2)';

trialDur = size(binData,2);
% betaMap = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
classLossTrain = zeros(trialDur/tBin,length(distStim));
classLossPredict = zeros(trialDur/tBin,length(distStim));
figure;axis square;
cMap = get(gca,'ColorOrder');

for iDist = 1:length(distStim)
    Cnt = 0;
    for x=1:tBin:trialDur
        Cnt = Cnt+1;
        
        correctInd = SessionData.StimType(perfTrials) == 3 & distFractions(perfTrials) == distStim(iDist) & SessionData.Rewarded(perfTrials); %multisensory trials at current distFrequency. Train on correct choice only
        trialCountTrain(Cnt,iDist) = sum(correctInd);
        SVMModel = fitcsvm(squeeze(mean(binData(:,x:x+tBin-1,correctInd),2))',leftInd(correctInd),'KernelFunction','linear','Standardize',true,'ClassNames',{'false','true'});
        %         betaMap(:,:,Cnt,iDist) = Widefield_DimMerge_v1(SVMModel.Beta,mask,'split'); %split pixels into original frame size
        
        CVSVMModel = crossval(SVMModel); %cross-validate mode;
        classLossTrain(Cnt,iDist) = kfoldLoss(CVSVMModel); %compute prediction error
        
        errorInd = SessionData.StimType(perfTrials) == 3 & distFractions(perfTrials) == distStim(iDist) & ~SessionData.Rewarded(perfTrials); %multisensory trials at current distFrequency. Test on error choice only
        trialCountPredict(Cnt,iDist) = sum(errorInd);
        classLossPredict(Cnt,iDist) = mean((ismember(predict(SVMModel,squeeze(mean(binData(:,x:x+tBin-1,errorInd),2))'),'true') - leftInd(errorInd)').^2); %try to predict outcome in error trials
        
        disp(['SVM map: Completed bin ' num2str(Cnt) '/' num2str(trialDur/tBin) '; DistFractions = ' num2str(distStim(iDist))]);
        %     SVMModel = fitclinear(squeeze(mean(binData(:,x:x+tBin-1,:),2)),iLeft,'ClassNames',{'false','true'},'ObserVationsIn','columns');
    end
    Est = [classLossTrain(1,iDist) time(round(length(time)*0.7)) (max(time)-min(time))/50 0.5];
    weights = ones(1,length(time)); 
%     weights(time>1 & time<=1.5) = 50;
    linearFitTrain{iDist} = fitnlm(time,classLossTrain(:,iDist),f,Est,'Weights',weights);
    
    %plot trained data
    subplot(2,length(distStim),iDist);
    plot(time,1-classLossTrain(:,iDist),'o','Markersize',10,'MarkerFaceColor','w','linewidth',2,'color',cMap(iDist,:));
    line(xx,1-predict(linearFitTrain{iDist},xx),'linewidth',2,'color',cMap(iDist,:)); % plot fit
    ylim([0.4 1]);xlim([-1 2]);title(['Dist. fraction: ' num2str(distStim(iDist))]);
    xlabel('time (s)');ylabel('Prediction error');axis square
    
    p = linearFitTrain{iDist}.Coefficients.Estimate';
    estAmps{iDist,1}(1) = p(2)-((p(3)*10)/2);    %M0 - onset of sigmoid rise
    estAmps{iDist,1}(2) = p(2)-((p(3)*2.2)/2);   %M25 - 25% of sigmoid rise
    estAmps{iDist,1}(3) = p(2);                  %M50 - inflection point
    estAmps{iDist,1}(4) = ((p(3)*2.2)/2)+p(2);   %M75 75% of sigmoid rise
    estAmps{iDist,1}(5) = ((p(3)*10)/2)+p(2);    %M100 - onset of sigmoid plateu
    estAmps{iDist,1}(estAmps{iDist,1} <= 0) = 0;  % amps can't be below 0
    line([estAmps{iDist,1}(3),estAmps{iDist,1}(3)],[-1 2],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
    
    %plot predicted data
    Est = [classLossPredict(1,iDist) time(round(length(time)*0.75)) (max(time)-min(time))/10 0.5];
    linearFitPredict{iDist} = fitnlm(time,classLossPredict(:,iDist),f,Est);

    subplot(2,length(distStim),iDist+length(distStim));
    plot(time,1-classLossPredict(:,iDist),'o','Markersize',10,'MarkerFaceColor','w','linewidth',2,'color',cMap(iDist,:));
    line(xx,1-predict(linearFitPredict{iDist},xx),'linewidth',2,'color',cMap(iDist,:)); % plot fit
    ylim([0.4 1]);xlim([-1 2]);title(['Dist. fraction: ' num2str(distStim(iDist))]);
    xlabel('time (s)');ylabel('Prediction error');axis square
    
    p = linearFitPredict{iDist}.Coefficients.Estimate';
    estAmps{iDist,2}(1) = p(2)-((p(3)*10)/2);    %M0 - onset of sigmoid rise
    estAmps{iDist,2}(2) = p(2)-((p(3)*2.2)/2);   %M25 - 25% of sigmoid rise
    estAmps{iDist,2}(3) = p(2);                  %M50 - inflection point
    estAmps{iDist,2}(4) = ((p(3)*2.2)/2)+p(2);   %M75 75% of sigmoid rise
    estAmps{iDist,2}(5) = ((p(3)*10)/2)+p(2);    %M100 - onset of sigmoid plateu
    estAmps{iDist,2}(estAmps{iDist,1} <= 0) = 0;  % amps can't be below 0
    line([estAmps{iDist,2}(3),estAmps{iDist,2}(3)],[-1 2],'linestyle','--','linewidth',1,'color',[0.5 0.5 0.5])
end

       
%% build choice classifier - perform random crossvalidation 100 times and average beta maps.
leftInd = SessionData.Assisted(perfTrials) & ...
    ((SessionData.CorrectSide(perfTrials) == 1 & SessionData.Rewarded(perfTrials)) | ... %animal went left and got rewarded
    (SessionData.CorrectSide(perfTrials) == 2 & SessionData.Punished(perfTrials))) ; %animal went left instead of right
% leftInd = SessionData.Assisted(perfTrials) & SessionData.CorrectSide(perfTrials) == 1; %leftward trials
tBin = 1;
crossValCnt = 10; %amount of repetitions for beta maps
holdOut = 0.5; %percent of data that should be used for crossvalidation

trialDur = size(binData,2);
betaMap = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
betaStd = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
classLoss = zeros(trialDur/tBin,length(distStim));
figure;axis square;hold on;
cMap = get(gca,'ColorOrder');

% for iDist = 1:length(distStim)
for iDist = length(distStim)
    Cnt = 0;
    for x=1:tBin:trialDur
        Cnt = Cnt+1;
        
        temp = zeros(size(binData,1),crossValCnt);
        for iModel = 1:crossValCnt
            
            trainSelect = rand(1,sum(perfTrials)) > holdOut; %randomly select subset of all trials for training
            detectInd = SessionData.StimType(perfTrials) == 3 & distFractions(perfTrials) >= distStim(iDist) & trainSelect; %multisensory trials at current distFrequency
            SVMModel = fitcsvm(squeeze(mean(binData(:,x:x+tBin-1,detectInd),2))',leftInd(detectInd),'KernelFunction','linear','Standardize',true);
            temp(:,iModel) = SVMModel.Beta; 
            testInd = SessionData.StimType(perfTrials) == 3 & distFractions(perfTrials) == distStim(iDist) & ~trainSelect; %use prediction subset to test current model
            modelLoss(iModel) = mean((leftInd(testInd)'-predict(SVMModel,squeeze(mean(binData(:,x:x+tBin-1,testInd),2))')).^2);
            
        end
        
        betaMapSingle(:,:,Cnt,iDist) = Widefield_DimMerge_v1(SVMModel.Beta,mask,'split'); %split pixels into original frame size
        betaMap(:,:,Cnt,iDist) = Widefield_DimMerge_v1(mean(temp,2),mask,'split'); %split pixels into original frame size
        betaStd(:,:,Cnt,iDist) = Widefield_DimMerge_v1(std(temp,0,2),mask,'split'); %split pixels into original frame size
        classLoss(Cnt,iDist) = mean(modelLoss); %compute prediction error
        
        disp(['SVM map: Completed bin ' num2str(Cnt) '/' num2str(trialDur/tBin) '; DistFractions = ' num2str(distStim(iDist))]);
        %     SVMModel = fitclinear(squeeze(mean(binData(:,x:x+tBin-1,:),2)),iLeft,'ClassNames',{'false','true'},'ObserVationsIn','columns');
    end
end

%% plot beta maps
% combinedBeta = abs(betaMap./betaStd); %normalize betaMap with standard deviaton

for iDist = 1:5
    figure
    subplot(1,3,1)
    imagesc(smooth2a(mean(abs(betaMap(:,:,:,iDist)),3),2));axis square;colormap jet
    title(['Distractor: ' num2str(distStim(iDist)) '; All frames'])
    subplot(1,3,2)
    imagesc(smooth2a(mean(abs(betaMap(:,:,80:120,iDist)),3),2));axis square;colormap jet
    title(['Distractor: ' num2str(distStim(iDist)) '; Frame 40-50'])
    subplot(1,3,3)
    imagesc(smooth2a(mean(abs(betaMap(:,:,120:160,iDist)),3),2));axis square;colormap jet
    title(['Distractor: ' num2str(distStim(iDist)) '; Frame 60-80'])
end


%% compute auc maps
%% build choice classifier - perform random crossvalidation 100 times and average beta maps.
leftInd = SessionData.Assisted(perfTrials) & ...
    ((SessionData.CorrectSide(perfTrials) == 1 & SessionData.Rewarded(perfTrials)) | ... %animal went left and got rewarded
    (SessionData.CorrectSide(perfTrials) == 2 & SessionData.Punished(perfTrials))) ; %animal went left instead of right
% leftInd = SessionData.Assisted(perfTrials) & SessionData.CorrectSide(perfTrials) == 1; %leftward trials
tBin = 1;
crossValCnt = 100; %amount of repetitions for beta maps
holdOut = 0.5; %percent of data that should be used for crossvalidation

trialDur = size(binData,2);
aucMap = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);
aucMapStd = zeros([size(mask),ceil(trialDur/tBin),length(distStim)]);

% for iDist = 1:length(distStim)
for iDist = [1 length(distStim)]
    Cnt = 0;
    for x=1:tBin:trialDur
        Cnt = Cnt+1;
        
        temp = zeros(size(binData,1),crossValCnt);
        for iModel = 1:crossValCnt
            trainSelect = rand(1,sum(perfTrials)) > holdOut; %randomly select subset of all trials for training
            detectInd = SessionData.StimType(perfTrials) == 3 & distFractions(perfTrials) == distStim(iDist) & trainSelect; %multisensory trials at current distFrequency
            temp(:,iModel) = colAUC((squeeze(mean(binData(:,x:x+tBin-1,detectInd),2))./4095)',leftInd(detectInd),'abs',false);
        if rem(iModel,10) == 0; disp(iModel);end
        end
        
        aucMap(:,:,Cnt,iDist) = Widefield_DimMerge_v1(mean(temp,2),mask,'split'); %split pixels into original frame size
        aucMapStd(:,:,Cnt,iDist) = Widefield_DimMerge_v1(std(temp,0,2),mask,'split'); %split pixels into original frame size
        
        disp(['AUC map: Completed bin ' num2str(Cnt) '/' num2str(trialDur/tBin) '; DistFractions = ' num2str(distStim(iDist))]);
        %     SVMModel = fitclinear(squeeze(mean(binData(:,x:x+tBin-1,:),2)),iLeft,'ClassNames',{'false','true'},'ObserVationsIn','columns');
    end
end
combinedAuc = abs(aucMap./aucMapStd); %normalize aucmap with standard deviaton


% trialDur = size(binData,2);
% tBin = 5;
% for iDist = 1:length(distStim)
%     Cnt = 0;
%     for x=1:tBin:trialDur
%         Cnt = Cnt+1;
%         detectInd = SessionData.StimType(perfTrials) == 3 & distFractions(perfTrials) == distStim(iDist); %multisensory trials at current distFrequency
%         aucMap(:,:,Cnt,iDist) = Widefield_DimMerge_v1(colAUC((squeeze(mean(binData(:,x:x+tBin-1,detectInd),2))./4095)',leftInd(detectInd),'abs',false)',mask,'split'); %compute area under the curve for ROC curves of all pixels and split auc values into original frame size
%         disp(['AUC map: Completed bin ' num2str(Cnt) '/' num2str(trialDur/tBin)  '; DistFractions = ' num2str(distStim(iDist))]);
%     end
% end
% aucMap = cat(3,Widefield_DimMerge_v1(squeeze(binData(:,1,1)),mask,'split'),aucMap);



%%
% figure
% for x=1:size(betaMap,3)
%     
%     %     subplot(1,2,1)
% %     imagesc(abs(squeeze(betaMap(:,:,x))));axis square; colormap jet; colorbar
%     imagesc(abs(smooth2a(squeeze(betaMap(:,:,x)),5,5)));axis square; colormap jet; colorbar
%     title(['Correct - ' num2str(x)]);
%     pause
% %     subplot(1,2,2)
% %     imagesc(squeeze(mean(binData(:,:,x,~iLeft),4)));axis square; colormap jet; caxis([0 0.01])
% %     title(['False - ' num2str(x)]);
% %     pause
%     
% end
% a

% 
% %% Look at single frames and check if there is a sustained response
% h1 = Widefield_ShowStack_v1(double(binData));
% 
% Widefield_SaveToAvi_v1(binData,[path Animal '\' Paradigm '\' Rec '\binData_5_Video'],5,'jet',[0 0.015])
% 
% if length(Frames) == 1
%     Check = false;
%     while ~Check
%         Wait = input('Enter specific frame nr. ','S');
%         if ~isnan(str2double(Wait))
%             Frames = str2double(Wait);
%             Check = true;
%         end
%     end
% end
% RejTrials = false;
% 
% %% Draw ROI through average and check if user is happy
% [aInd,mask,h2] = FindROI_v4(binData,Frames(Frames>=stimOn),2,rThresh,smth,1,'k',cRange); %find ROI
% 
% Check = false;
% while ~Check
%     Wait = input('Happy with ROI? (Y/N)  ','S');
%     if strcmpi(Wait,'y')
%         Check = true;
%     elseif strcmpi(Wait,'n')
%         close(h2);
%         [aInd,mask,h2] = FindROI_v4(binData,Frames,2,rThresh,smth,1,'k',cRange); %find ROI
%     elseif  ~isempty(str2num(Wait))
%         close(h2);
%         Wait = str2num(Wait);
%         if length(Wait) == 1
%             disp(['Changed threshold to '  num2str(Wait)]);
%             rThresh = Wait; %treshold for ROI
%         elseif length(Wait) == 2
%             disp(['Changed caxis to ' num2str(Wait)]);
%             cRange = Wait; %treshold for ROI
%         end
%         [aInd,mask,h2] = FindROI_v4(binData,Frames,2,rThresh,smth,1,'k',cRange); %find ROI
%     end
% end
% 
% %% If result is insufficient, try to reject bad trials
% RedoRejection = true;
% tInd = true(1,size(binData,4)); %use all trials
% 
% if RejTrials
%     if exist([path Animal '\' Paradigm '\' Rec '\RejTrials.mat'],'file')
%         Check = false;
%         while ~Check
%             a = input('Reuse old rejection vector? (Y/N) ','S');
%             if strcmpi(a,'y')
%                 load([path Animal '\' Paradigm '\' Rec '\RejTrials.mat'])
%                 Check = true; RedoRejection = false;
%             elseif strcmpi(a,'n')
%                 Check = true; RedoRejection = true;
%             end
%         end
%     end
%     
%     if RedoRejection
%         close(h2);
%         h2 = figure;colormap('jet')
%         for x=1:size(binData,4)
%             imagesc(squeeze(mean(binData(:,:,(Frames),x),3)));axis square;caxis([-0.015  0.015])
%             Check = false;
%             while ~Check
%                 UseTrial = input(['Use trial ' int2str(x) ' ? '],'S');
%                 if strcmpi(UseTrial,'Y')
%                     UseTrial = true;
%                     Check = true;
%                 elseif strcmpi(UseTrial,'N')
%                     UseTrial = false;
%                     Check = true;
%                 end
%             end
%             tInd(x) = UseTrial;
%         end
%     end
%     close(h2);
%     [aInd,mask,~] = FindROI_v4(binData(:,:,:,tInd),Frames,2,rThresh,smth,1,'k'); %find ROI
%     save([path Animal '\' Paradigm '\' Rec '\RejTrials.mat'],'tInd')
% end
% 
% %% get index for selected ROI and compute activity traces
% xInd = [];
% yInd = [];
% xFind = find(sum(mask,1)>0);
% 
% for x = 1:length(xFind)
%     yFind = find(mask(:,xFind(x)));
%     for y = 1:length(yFind)
%         xInd = [xInd xFind(x)];
%         yInd = [yInd yFind(y)];
%     end
% end
% 
% Cnt = 1;
% for x = 1:length(xInd)
%     roiTrace(Cnt,:,:) = binData(yInd(x),xInd(x),:,:); %activity trace in selected ROI
%     Cnt = Cnt+1;
% end
% wholeTrace = squeeze(nanmean(nanmean(nanmean(binData,1),2),4)); %activity trace of all frames
%     
% 
% %% load vessel image and add ROI to it
% % figure;
% % imagesc(snap); axis square; hold on; colormap gray
% % plot(aInd(1,:)*binSize,aInd(2,:)*binSize,'c','linewidth',2)
% % title([Animal ' - ' Rec])
% 
% figure;
% subplot(2,2,1)
% imagesc(snap); axis square; hold on; colormap(gca,'gray')
% plot(aInd(1,:),aInd(2,:),'k','linewidth',2)
% title([Animal ' - ' Rec])
% 
% subplot(2,2,2)
% title([Animal ' - ' Rec])
% imagesc(smooth2a(nanmean(nanmean(binData(:,:,Frames,tInd),3),4),smth,smth));axis square;colormap(gca,'jet')
% hold on;plot(aInd(1,:),aInd(2,:),'k','linewidth',2)
% title([Animal ' - ' Rec])
% 
% subplot(2,2,3)
% plot(wholeTrace);axis square;title('Whole image activity')
% xlabel('# Frame');ylabel('delta R / R')
% 
% subplot(2,2,4)
% plot((squeeze(nanmean(nanmean(roiTrace,1),3))));axis square;title('ROI activity')
% xlabel('# Frame')
% 
% %%
% close(h1)