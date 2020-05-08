function Widefield_testLinearModel(cPath,Animal,Rec)

%general variables
Paradigm = 'SpatialDisc';
cPath = [cPath Animal filesep Paradigm filesep Rec filesep]; %Widefield data path
% sPath = ['U:\space_managed_data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
% sPath = ['/sonas-hs/churchland/hpc/home/space_managed_data/BpodImager/Animals/' Animal filesep Paradigm filesep Rec filesep]; %server data path. not used on hpc.
ridgeFolds = 10; %number of folds for cross-validation
dims = 200; %number of dimensions, used from Vc

%% load imaging data and design matrix
load([cPath filesep 'regData.mat']);
load([cPath filesep 'ridgeTest.mat']);
load([cPath filesep 'Vc.mat']);
load([cPath filesep 'mask.mat']);

[~,frames,~] = size(Vc);
Vc = reshape(Vc(1:dims,:),dims,[]);
Vc(:,trialIdx) = [];
U = arrayShrink(U(:,:,1:dims),mask);

[~, R] = qr(fullR./sqrt(sum(fullR.^2)),0);
if ~fullColRank(R)
    error('Design matrix is not orthogonal enough. Major danger !')
end

%% cycle through regressor sets and compare vs reduced model
Vc = reshape(Vc,dims,[]);

% allocate larger data arrays
meanPredTrial = zeros(length(unique(recIdx))+1,frames);
semPredTrial = zeros(length(unique(recIdx))+1,frames);
corrMap = zeros(length(unique(recIdx))+1,sum(~mask(:)));
maxCorrMap = zeros(length(unique(recIdx))+1,sum(~mask(:)));

Cnt = 0;
for iRegs = [0 unique(recIdx)]
    
    Cnt = Cnt+1;
    fprintf('Current regressor is %d out of %d\n', Cnt,length(unique(recIdx)+1));
    
    if iRegs == 0
        cIdx = true(1,size(fullR,2)); %use all regressors for full model
    else
        cIdx = recIdx(~idx) ~= iRegs; %index for reduced model.
    end    
    
    fakeR = fullR; %create shuffled comparison dataset for current regressor

    [meanRsq, trialRsq, cMap, mcMap] = Widefield_crossValModel(Vc,fullR,U,cIdx,ridgeFolds,frames);

    meanPredAll(Cnt) = mean(meanRsq);
    semPredAll(Cnt) = sem(meanRsq);
    meanPredTrial(Cnt,:) = mean(trialRsq);
    semPredTrial(Cnt,:) = sem(trialRsq);
    corrMap(Cnt,:) = cMap;
    maxCorrMap(Cnt,:) = mcMap;
    
    fakeR = fullR; %create shuffled comparison dataset for current regressor
    fakeR(:,~cIdx) = fullR(randperm(size(fullR,1)),~cIdx); %shuffle rows of current regressors around
    
    fakeRsq = Widefield_crossValModel(Vc,fullR,U,cIdx,ridgeFolds,frames);
    meanPredFake(Cnt) = mean(fakeRsq);
    semPredFake(Cnt) = sem(fakeRsq);
    
    [rankP(Cnt), rankH(Cnt)] = signrank(meanRsq, meanPredFake(Cnt),'tail', 'left'); %check for significant reduction vs full model
    fprintf('Finished. Current reg: %d ; pVal is %d\n', Cnt,rankP(Cnt));

end

%% save results
recLabels = {'full' recLabels{:}};
save([cPath filesep 'modelCompare.mat'],'meanPredAll','semPredAll','meanPredTrial','semPredTrial', 'corrMap', 'maxCorrMap', 'meanPredFake', 'semPredFake', 'rankP', 'rankH','recLabels');

%%

[dMean,b] = sort(meanPredAll,'descend');
dSem = semPredAll(b);

figure
ax = stairs(dMean,'o','lineWidth',2,'MarkerSize',10);
% vline(find(h(b)))
hline(mean(dMean(1)))
set(ax.Parent,'xTick',1:iRegs)
set(ax.Parent,'xTickLabel',recLabels(b))
