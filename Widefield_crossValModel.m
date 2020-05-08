function [trialRsq, totalRsq, timeRsq, corrMap, corrMovie, Vm] = Widefield_crossValModel(Vc,fullR,U,ridgeFolds,frames, useGPU, makeMovie,randSampling)

if ~exist('useGPU','var')
    useGPU = false; %flag to use GPU
end

if ~exist('makeMovie','var')
    makeMovie = false; %flag to produce correlation movie
end

if ~exist('randSampling','var')
    randSampling = false; %flag to use a random instead of linear index for cross-validation
end

%% run cross-validation
if randSampling
    rng default %reset randum number generator
    randIdx = randperm(size(Vc,2)); %generate randum number index if required
end

Vm = zeros(size(Vc),'single');
foldCnt = floor(size(Vc,2) / ridgeFolds);

tic
for iFolds = 1:ridgeFolds
    tic
    dataIdx = true(1,size(Vc,2));
    
    if ridgeFolds > 1
        
        if randSampling
            dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
        else
            dataIdx(((iFolds - 1)*foldCnt) + (1:foldCnt)) = false; %index for training data
        end
        
        Y = Vc(:,dataIdx)';
        Y = bsxfun(@minus,Y,mean(Y));
        
        X = fullR(dataIdx,:);
        Xm = fullR(~dataIdx,:);
        
        X = bsxfun(@minus,X,mean(X));
        Xm = bsxfun(@minus,Xm,mean(X));
        
        [~, betas] = ridgeMML(Y, X); %get beta weights for current model
        Vm(:,~dataIdx) = (Xm * betas)'; %predict remaining data
        
        if rem(iFolds,ridgeFolds/5) == 0
            fprintf(1, 'Current fold is %d out of %d\n', iFolds,ridgeFolds);
            toc
        end
    else
        
        Vc = bsxfun(@minus,Vc,mean(Vc));
        fullR = bsxfun(@minus,fullR,mean(fullR));
        
        [~, betas] = ridgeMML(Vc', fullR); %get beta weights for current model
        Vm = (fullR * betas)'; %predict remaining data
        disp('Ridgefold is <= 1, fit to complete dataset');
    end
end

Vc = reshape(Vc,size(Vc,1),frames,[]);
Vm = reshape(Vm,size(Vm,1),frames,[]);

cError = (Vc-Vm).^2;

totalRsq = 1 - (mean(cError(:)) ./ var(Vc(:))); %total mean R^2
trialRsq = squeeze(1 - (mean(mean(cError),2) ./ mean(var(Vc,[],2),1))); %mean R^2 per trial
timeRsq = 1 - (squeeze(mean(mean(cError,3),1) ./ mean(var(Vc,[],3),1))); %mean R^2 per frame / trial

%% compute correlation between data and prediction over for all frames
if useGPU
    Vc = gpuArray(Vc);
    Vm = gpuArray(Vm);
    U = gpuArray(U);
end

tic
Vc = reshape(Vc,size(Vc,1),[]);
Vm = reshape(Vm,size(Vm,1),[]);
covVc = cov(Vc');  % S x S
covVm = cov(Vm');  % S x S
cCovV = bsxfun(@minus, Vm, mean(Vm,2)) * Vc' / (size(Vc, 2) - 1);  % S x S
covP = sum((U * cCovV) .* U, 2)';  % 1 x P
varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
corrMap = gather((covP ./ stdPxPy)');
toc

if makeMovie
    for iFrames = 1:frames
        idx = iFrames:frames:size(Vc,2); %index for frame with lowest prediction error in each trial
        tic
        covVc = cov(Vc(:,idx)');  % S x S
        covVm = cov(Vm(:,idx)');  % S x S
        cCovV = bsxfun(@minus, Vm(:,idx), mean(Vm(:,idx),2)) * Vc(:,idx)' / (length(idx) - 1);  % S x S
        covP = sum((U * cCovV) .* U, 2)';  % 1 x P
        varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
        varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
        stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
        corrMovie(:,iFrames) = (covP ./ stdPxPy)';
        
        if rem(iFrames,frames/20) == 0
            fprintf(1, 'Current frame is %d out of %d\n', iFrames,frames);
            toc
        end
    end
    corrMovie = gather(corrMovie);
else
    corrMovie = corrMap;
end
fprintf('Run finished. RMSE: %f\n', mean(totalRsq));

end