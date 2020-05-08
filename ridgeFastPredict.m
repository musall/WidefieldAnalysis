function [meanRMSE, ridgeRMSE] = ridgeFastPredict(Vc, fullR, ridgeFolds, ridgeVal)
% perform ridge regression with xRidgeFolds cross-validation.
% Returns RMSE for optimization to identify optimal ridge-penalty.
    
foldCnt = floor(size(Vc,2) / ridgeFolds);
for iFolds = 1:ridgeFolds
    cIdx = true(1,size(Vc,2));
    cIdx(((iFolds - 1)*foldCnt) + (1:foldCnt)) = false;
    [~, ridgeRMSE(iFolds)] = ridgeFast(Vc', fullR, ridgeVal, cIdx);
end
meanRMSE = mean(ridgeRMSE);
fprintf('Run finished. RidgeVal: %f, RMSE: %f\n', ridgeVal, meanRMSE);