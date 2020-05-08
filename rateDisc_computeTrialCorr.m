function [trialVar, timeVar] = rateDisc_computeTrialCorr(U, Vc, Vm, mask, allenMask, opts)
% function to compute explained variance for a given model Vm and real data Vm. 
% Assumes that U is based on mask and realigns the output to the allen using opts. 

% compute trial-by-trial correlations
timeAvg = squeeze(mean(Vc,2));
corrMat = rateDisc_modelCorr(timeAvg,squeeze(mean(Vm,2)),U) .^2; %R2 map for trial-by-trial variability
corrMat = arrayShrink(corrMat,mask,'split');
corrMat = alignAllenTransIm(corrMat,opts.transParams);
trialVar = arrayShrink(corrMat(:,1:size(allenMask,2)),allenMask, 'merge');

% compute trial-averaged correlations
trialAvg = mean(Vc,3);
corrMat = rateDisc_modelCorr(trialAvg,mean(Vm,3),U) .^2; %R2 map for trial-by-trial variability
corrMat = arrayShrink(corrMat,mask,'split');
corrMat = alignAllenTransIm(corrMat,opts.transParams);
timeVar = arrayShrink(corrMat(:,1:size(allenMask,2)),allenMask, 'merge');