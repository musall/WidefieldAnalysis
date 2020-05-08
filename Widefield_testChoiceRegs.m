function Widefield_testChoiceRegs(fPath, Animal, Rec, cReg)
% Code to compute predictive power in different choice regressors

if ~strcmpi(fPath(end),filesep)
    fPath = [fPath filesep];
end

fPath = [fPath Animal filesep 'SpatialDisc' filesep Rec filesep]; %Widefield data path

%% load some data
[~, motorLabels, sensorLabels, cogLabels] = delayDecRecordings;
load([fPath 'mask.mat'], 'mask')
load([fPath 'Vc.mat'], 'U')
load([fPath 'orthChoicedimBeta.mat'])
load([fPath 'orthChoiceregData.mat'])
load([fPath 'ninterpVc.mat'], 'Vc', 'frames')
U = arrayShrink(U, mask);

%indices for trial segments. Baseline and handle are based on time regressors, remaining segments are based on stimulus regressors.
segIdx = {1:18 55:72 1:18 34:51 57:74 81:98}; %set index so different segments have same #frames

ridgeFolds = 10;    %folds for cross-validation when assessing predicted variance
rng default %reset randum number generator
randIdx = randperm(size(Vc,2)); %generate randum number index if required

%% test predictive power of choice regressors through cross-validation
cIdx = ismember(recIdx, find(ismember(recLabels,cReg)));
fakeR = fullR; %copy  design matrix to shuffle up some regressor set

%shuffle selected regressors
for iCol = find(cIdx)
    fakeR(:,iCol) = fullR(randperm(size(fullR,1)),iCol);
end

%% run cross-validation
Vm = zeros(size(Vc),'single');
foldCnt = floor(size(Vc,2) / ridgeFolds);

tic
for iFolds = 1:ridgeFolds
    tic
    dataIdx = true(1,size(Vc,2));
    
    if ridgeFolds > 1
        
        dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
        
        [~, betas] = ridgeMML(Vc(:,dataIdx)', fakeR(dataIdx,:), true, ridgeVals); %get beta weights for current model. Re-center and supply ridge penalty.
        Vm(:,~dataIdx) = (fakeR(~dataIdx,:) * betas)'; %predict remaining data
        
        if rem(iFolds,ridgeFolds/5) == 0
            fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
            toc
        end
    else
        
        [~, betas] = ridgeMML(Vc', fakeR, true, ridgeVals); %get beta weights for current model. Re-center and supply ridge penalty.
        Vm = (fakeR * betas)'; %predict remaining data
        disp('Ridgefold is <= 1, fit to complete dataset');
    end
end

%% compute correlation between data and prediction over for all frames
if ispc
    Vc = gpuArray(Vc);
    Vm = gpuArray(Vm);
    U = gpuArray(U);
end

Vc = reshape(Vc,size(Vc,1),[]);
Vm = reshape(Vm,size(Vm,1),[]);
covVc = cov(Vc');  % S x S
covVm = cov(Vm');  % S x S
cCovV = bsxfun(@minus, Vm, mean(Vm,2)) * Vc' / (size(Vc, 2) - 1);  % S x S
covP = sum((U * cCovV) .* U, 2)';  % 1 x P
varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
cMap = gather((covP ./ stdPxPy)');

%% compute correlation between data and prediction over specific trial segments
stimRegs = sum(fullR(:,ismember(recIdx,find(ismember(recLabels, sensorLabels)))),2); %index for stimulus onset in all trials
stimRegs = find([0;diff(stimRegs)] > 0.5) - 1;

segMovie = zeros(size(U,1),length(segIdx), 'single');
for iSegs = 1:length(segIdx)
    if iSegs <= 2 %baseline and handle grab, use trial onset times
        cIdx = repmat((0:frames:size(Vc,2)-1)',1,length(segIdx{iSegs}));
    else %stimulus and later segments, use stimulus onset times
        cIdx = repmat(stimRegs,1,length(segIdx{iSegs}));
    end
    cIdx = bsxfun(@plus,cIdx,segIdx{iSegs});
    cIdx = cIdx(:);
    cIdx(cIdx > size(Vc,2)) = [];
    
    covVc = cov(Vc(:,cIdx)');  % S x S
    covVm = cov(Vm(:,cIdx)');  % S x S
    cCovV = bsxfun(@minus, Vm(:,cIdx), mean(Vm(:,cIdx),2)) * Vc(:,cIdx)' / (length(cIdx) - 1);  % S x S
    covP = sum((U * cCovV) .* U, 2)';  % 1 x P
    varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
    varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
    stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
    segMovie(:,iSegs) = gather((covP ./ stdPxPy))';
end

%% compute correlation between data and prediction for each frame in all trials
cMovie = zeros(size(U,1),frames, 'single');
for iFrames = 1:frames
    
    frameIdx = iFrames:frames:size(Vc,2); %index for the same frame in each trial
    tic
    cData = bsxfun(@minus, Vc(:,frameIdx), mean(Vc(:,frameIdx),2));
    cModel = bsxfun(@minus, Vm(:,frameIdx), mean(Vm(:,frameIdx),2));
    covVc = cov(cData');  % S x S
    covVm = cov(cModel');  % S x S
    cCovV = cModel * cData' / (length(frameIdx) - 1);  % S x S
    covP = sum((U * cCovV) .* U, 2)';  % 1 x P
    varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
    varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
    stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
    cMovie(:,iFrames) = gather(covP ./ stdPxPy)';
    clear cData cModel
    
    if rem(iFrames,round(frames/4)) == 0
        fprintf(1, 'Current frame is %d out of %d\n', iFrames,frames);
        toc
    end
end
fprintf('Run finished. RMSE: %f\n', median(cMovie(:).^2));

%% save results
if ~exist([fPath 'predVariance' filesep 'newShCurrent'],'dir')
    mkdir([fPath 'predVariance' filesep 'newShCurrent']);
end
save([fPath 'predVariance' filesep 'newShCurrent' filesep cReg 'corr.mat'], 'cMap', 'cMovie', 'segMovie', 'recLabels', 'cReg', '-v7.3');
end