function delayDec_batchNoAnalog(cMode)

if ~exist('cMode')
    cMode = 'show';
end

% cPath = 'Y:\data\BpodImager\Animals\';
% cPath = 'X:\smusall\BpodImager\Animals\';
iMod = 3;
if iMod == 1
    % GFP data
    [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
    dataOverview = dataOverview(ismember(dataOverview(:,4), 'GFP'), :);
elseif iMod == 2
    % camK data
    [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
    dataOverview = dataOverview(ismember(dataOverview(:,4), 'gcamp'), :);
    cPath = 'U:\smusall\BpodImager\Animals\'; %Widefield data path on grid server
elseif iMod == 3
    % ai93 data
    [dataOverview, motorLabels, ~, ~, ~, ~, ~, cPath] = delayDecRecordings;
end
load('allenDorsalMapSM.mat','dorsalMaps')
allenMask = dorsalMaps.allenMask; %mask that is used for all datasets
cPath = 'U:\smusall\BpodImager\Animals\';
animals = dataOverview(:,1);
recs = dataOverview(:,3);
opMotorLabels = {'lLick' 'rLick' 'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel'}; %operant motor regressors
spontLabels = motorLabels(~ismember(motorLabels,opMotorLabels));
ridgeFolds = 10;

corrMaps = NaN(sum(~allenMask(:)), 7, length(animals), 'single');
for iAnimals = 1 : length(animals)
    tic
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep];
    disp(fPath); %current data path
    
    %%
    if strcmpi(cMode, 'reload')
        %% load design matrix and reject analog regressors
        load([fPath 'orgregData.mat'], 'fullR', 'recIdx', 'idx', 'recLabels')
        load([fPath 'Vc.mat'], 'U')
        load([fPath 'interpVc.mat'], 'Vc', 'frames')
        recIdx = recIdx(~idx); clear idx
        
        rejIdx = false(1, size(fullR,2));
        for iRegs = 1 : size(fullR,2)
            if length(unique(fullR(:,iRegs))) > 2
                rejIdx(iRegs) = true;
            end
        end
        fullR(:,rejIdx) = [];
        recIdx(rejIdx) = [];
        disp(sum(rejIdx));
        
        %% orthogonalize some regressors for clarity
        % orthogonalize spontaneous from operant movement regressors
        lInd = ismember(recIdx, find(ismember(recLabels,{'lLick', 'rLick'})));
        hInd = ismember(recIdx, find(ismember(recLabels,{'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel'})));
        pInd = ismember(recIdx, find(ismember(recLabels,{'fastPupil', 'slowPupil'})));
        wInd = ismember(recIdx, find(ismember(recLabels,'whisk')));
        nInd = ismember(recIdx, find(ismember(recLabels,'nose')));
        piInd = ismember(recIdx, find(ismember(recLabels,'piezo')));
        fInd = ismember(recIdx, find(ismember(recLabels,'face')));
        
        smallR = [fullR(:,lInd) fullR(:,hInd) fullR(:,pInd) fullR(:,wInd) fullR(:,nInd) fullR(:,piInd) fullR(:,fInd)];
        [Q, ~] = qr(smallR,0); clear smallR %orthogonalize spont. from operant movement
        
        % replace original with orthogonalized regressors (only for spont. movements)
        fullR(:,pInd) = Q(:,sum(lInd | hInd) + 1 : sum(lInd | hInd | pInd)); %pupil
        fullR(:,wInd) = Q(:,sum(lInd | hInd | pInd) + 1 : sum(lInd | hInd | pInd | wInd)); %whisk
        fullR(:,nInd) = Q(:,sum(lInd | hInd | pInd | wInd) + 1 : sum(lInd | hInd | pInd | wInd | nInd)); %nose
        fullR(:,piInd) = Q(:,sum(lInd | hInd | pInd | wInd | nInd) + 1 : sum(lInd | hInd | pInd | wInd | nInd | piInd)); %piezo
        fullR(:,fInd) = Q(:,sum(lInd | hInd | pInd | wInd | nInd | piInd) + 1 : sum(lInd | hInd | pInd | wInd | nInd | piInd | fInd)); %face
        
        %%
        [fVm, fBeta, fR, fIdx, fRidge, frecLabels, fcMap, fcMovie] = crossValModel(recLabels);
        save([fPath 'filterVm.mat'], 'fVm', 'frames'); %save predicted data based on full filter-only model
        save([fPath 'filterBeta.mat'], 'fBeta', 'fRidge');
        save([fPath 'filterregData.mat'], 'fR','fIdx', 'frecLabels', 'fcMap', 'fcMovie', '-v7.3');
        
        [Vmotor, motorBeta, motorR, motorIdx, motorRidge, motorLabels, motorMap, motorMovie] = crossValModel(motorLabels);
        save([fPath 'filterVmotor.mat'], 'Vmotor', 'frames'); %save predicted data based on motor model
        save([fPath 'filterMotorBeta.mat'], 'motorBeta', 'motorRidge');
        save([fPath 'filterMotorregData.mat'], 'motorR', 'motorIdx', 'motorLabels', 'motorMap', 'motorMovie', '-v7.3');

        [Vtask, taskBeta, taskR, taskIdx, taskRidge, taskLabels, taskMap, taskMovie] = crossValModel(recLabels(~ismember(recLabels,motorLabels)));
        save([fPath 'filterVtask.mat'], 'Vtask', 'frames'); %save predicted data based on task model
        save([fPath 'filterTaskBeta.mat'], 'taskBeta', 'taskRidge');
        save([fPath 'filterTaskregData.mat'], 'taskR', 'taskIdx', 'taskLabels', 'taskMap', 'taskMovie', '-v7.3');

        [VtaskOpMotor, taskOpMotorBeta, taskOpMotorR, taskOpMotorIdx, taskOpMotorRidge, taskOpMotorLabels, taskOpMotorMap, taskOpMotorMovie] = crossValModel([recLabels(~ismember(recLabels,motorLabels)) opMotorLabels]);
        save([fPath 'filterVtaskOpMotor.mat'], 'VtaskOpMotor', 'frames'); %save predicted data based on taskOpMotor model
        save([fPath 'filterTaskOpMotorBeta.mat'], 'taskOpMotorBeta', 'taskOpMotorRidge');
        save([fPath 'filterTaskOpMotorregData.mat'], 'taskOpMotorR','taskOpMotorIdx', 'taskOpMotorLabels', 'taskOpMotorMap', 'taskOpMotorMovie', '-v7.3');
        
        [VtaskSpMotor, taskSpMotorBeta, taskSpMotorR, taskSpMotorIdx, taskSpMotorRidge, taskSpMotorLabels, taskSpMotorMap, taskSpMotorMovie] = crossValModel(recLabels(~ismember(recLabels,opMotorLabels)));
        save([fPath 'filterTaskSpMotor.mat'], 'VtaskSpMotor', 'frames'); %save predicted data based on taskSpMotor model
        save([fPath 'filterTaskSpMotorBeta.mat'], 'taskSpMotorBeta', 'taskSpMotorRidge');
        save([fPath 'filterTaskSpMotorregData.mat'], 'taskSpMotorR','taskSpMotorIdx', 'taskSpMotorLabels', 'taskSpMotorMap', 'taskSpMotorMovie', '-v7.3');

        [VspontMotor, spontMotorBeta, spontMotorR, spontMotorIdx, spontMotorRidge, spontMotorLabels, spontMotorMap, spontMotorMovie] = crossValModel(motorLabels(~ismember(motorLabels,opMotorLabels)));
        save([fPath 'filterVspontMotor.mat'], 'VspontMotor', 'frames'); %save predicted data based on spontaneous motor model
        save([fPath 'filterSpontMotorBeta.mat'], 'spontMotorBeta', 'spontMotorRidge');
        save([fPath 'filterSpontMotorregData.mat'], 'spontMotorR', 'spontMotorIdx', 'spontMotorLabels', 'spontMotorMap', 'spontMotorMovie','-v7.3');
        
        [VopMotor, opMotorBeta, opMotorR, opMotorIdx, opMotorRidge, opMotorLabels, opMotorMap, opMotorMovie] = crossValModel(opMotorLabels);
        save([fPath 'filterVopMotor.mat'], 'VopMotor', 'frames'); %save predicted data based on operant motor model
        save([fPath 'filterOpMotorBeta.mat'], 'opMotorBeta', 'opMotorRidge');
        save([fPath 'filterOpMotorregData.mat'], 'opMotorR', 'opMotorIdx', 'opMotorLabels', 'opMotorMap', 'opMotorMovie', '-v7.3');
        
    elseif strcmpi(cMode, 'show')
        
        load([fPath 'filterregData.mat'], 'fcMap');
        load([fPath 'filterMotorregData.mat'], 'motorMap');
        load([fPath 'filterTaskregData.mat'], 'taskMap');
        load([fPath 'filterTaskOpMotorregData.mat'], 'taskOpMotorMap');
        load([fPath 'filterTaskSpMotorregData.mat'], 'taskSpMotorMap');
        load([fPath 'filterSpontMotorregData.mat'], 'spontMotorMap');
        load([fPath 'filterOpMotorregData.mat'], 'opMotorMap');

    end
    
    load([fPath 'mask.mat'], 'mask');
    load([fPath 'opts2.mat'], 'opts');
    
    cMap = arrayShrink(fcMap.^2, mask, 'split');
    cMap = alignAllenTransIm(cMap,opts.transParams);
    corrMaps(:, 1, iAnimals) = arrayShrink(cMap(1:size(allenMask,1),1:size(allenMask,2)), allenMask, 'merge');
    
    cMap = arrayShrink(taskSpMotorMap.^2, mask, 'split');
    cMap = alignAllenTransIm(cMap,opts.transParams);
    corrMaps(:, 2, iAnimals) = arrayShrink(cMap(1:size(allenMask,1),1:size(allenMask,2)), allenMask, 'merge');
    
    cMap = arrayShrink(taskOpMotorMap.^2, mask, 'split');
    cMap = alignAllenTransIm(cMap,opts.transParams);
    corrMaps(:, 3, iAnimals) = arrayShrink(cMap(1:size(allenMask,1),1:size(allenMask,2)), allenMask, 'merge');
        
    cMap = arrayShrink(motorMap.^2, mask, 'split');
    cMap = alignAllenTransIm(cMap,opts.transParams);
    corrMaps(:, 4, iAnimals) = arrayShrink(cMap(1:size(allenMask,1),1:size(allenMask,2)), allenMask, 'merge');
            
    cMap = arrayShrink(taskMap.^2, mask, 'split');
    cMap = alignAllenTransIm(cMap,opts.transParams);
    corrMaps(:, 5, iAnimals) = arrayShrink(cMap(1:size(allenMask,1),1:size(allenMask,2)), allenMask, 'merge');
                
    cMap = arrayShrink(spontMotorMap.^2, mask, 'split');
    cMap = alignAllenTransIm(cMap,opts.transParams);
    corrMaps(:, 6, iAnimals) = arrayShrink(cMap(1:size(allenMask,1),1:size(allenMask,2)), allenMask, 'merge');
    
    cMap = arrayShrink(opMotorMap.^2, mask, 'split');
    cMap = alignAllenTransIm(cMap,opts.transParams);
    corrMaps(:, 7, iAnimals) = arrayShrink(cMap(1:size(allenMask,1),1:size(allenMask,2)), allenMask, 'merge');
    
end

%% make overview figure for overall task- and motor explained variance
figure;
cRange = 0.1;
figure;
subplot(2,2,1);
fullMap = arrayShrink(nanmean(corrMaps(:,1,:),3),allenMask,'split');
mapImg = imshow(fullMap,[0 0.5]);axis image;
set(mapImg,'AlphaData',~allenMask); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256));colorbar; title('Full model')
    
subplot(2,2,2);
mapImg = imshow(arrayShrink(nanmean(corrMaps(:,1,:) - corrMaps(:,3,:),3),allenMask,'split'),[0 cRange]);axis image;
set(mapImg,'AlphaData',~allenMask); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256));colorbar; title('Uninstructed movements')

subplot(2,2,3);
mapImg = imshow(arrayShrink(nanmean(corrMaps(:,1,:) - corrMaps(:,4,:),3),allenMask,'split'),[0 cRange]);axis image;
set(mapImg,'AlphaData',~allenMask); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256));colorbar; title('Task model')

subplot(2,2,4);
mapImg = imshow(arrayShrink(nanmean(corrMaps(:,1,:) - corrMaps(:,2,:),3),allenMask,'split'),[0 cRange]);axis image;
set(mapImg,'AlphaData',~allenMask); %make NaNs transparent.
colormap(mapImg.Parent,inferno(256)); colorbar; title('Instructed movements')

disp(['Full model performance: ' num2str(nanmean(fullMap(:))*100) ' ± ' num2str(sem(fullMap(:))*100) ' %']);
    
%% make figure - R2 bars
cData = [];
cData(1,:,:) = squeeze(nanmean(corrMaps(:,[5 6 7],:)))';
cData(2,:,:) = squeeze(nanmean(corrMaps(:,1,:)) - nanmean(corrMaps(:,[4 3 2],:)))';

figure
[ax, idxGroup] = regressorBoxPlot(squeeze(cData(1,:,:)), {'Task' 'spontMove' 'opMove'}, 5, subplot(2,2,1:4), [0 1 0], {[1 3 2]},0.6,[-0.2 0.65]);
ax = regressorBoxPlot(squeeze(-cData(2,:,:)), {'Task' 'spontMove' 'opMove'}, 5, ax, [25 111 61]/255, idxGroup,0.6,[-0.2 0.65]);
title('Regg groups: All (top) / Unique (bottom)'); xlim([0 30]); %make this a fixed value to ensure bars have always the same width
ylabel('total cv R^2');
hline(nanmean(nanmean(corrMaps(:,1,:))),'k--')


%%
function [Vm, cBeta, cR, subIdx, cRidge, cLabels, cMap, cMovie] =  crossValModel(cLabels)

cIdx = ismember(recIdx, find(ismember(recLabels,cLabels))); %get index for task regressors
cLabels = recLabels(sort(find(ismember(recLabels,cLabels)))); %make sure motorLabels is in the right order

%create new regressor index that matches motor labels
subIdx = recIdx;
subIdx = subIdx(cIdx);
temp = unique(subIdx);
for x = 1 : length(temp)
    subIdx(subIdx == temp(x)) = x;
end
cR = fullR(:,cIdx);

Vm = zeros(size(Vc),'single'); %pre-allocate motor-reconstructed V
randIdx = randperm(size(Vc,2)); %generate randum number index
foldCnt = floor(size(Vc,2) / ridgeFolds);
cBeta = cell(1,ridgeFolds);

for iFolds = 1:ridgeFolds
    dataIdx = true(1,size(Vc,2));
    
    if ridgeFolds > 1
        dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
        if iFolds == 1
            [cRidge, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', cR(dataIdx,:), true); %get beta weights and ridge penalty for task only model
        else
            [~, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', cR(dataIdx,:), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
        end
        Vm(:,~dataIdx) = (cR(~dataIdx,:) * cBeta{iFolds})'; %predict remaining data
        
        if rem(iFolds,ridgeFolds/5) == 0
            fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
        end
    else
        [cRidge, cBeta{iFolds}] = ridgeMML(Vc', cR, true); %get beta weights for task-only model.
        Vm = (cR * cBeta{iFolds})'; %predict remaining data
        disp('Ridgefold is <= 1, fit to complete dataset instead');
    end
end

% computed all predicted variance
Vc = reshape(Vc,size(Vc,1),[]);
Vm = reshape(Vm,size(Vm,1),[]);
if length(size(U)) == 3
    U = arrayShrink(U, squeeze(isnan(U(:,:,1))));
end
covVc = cov(Vc');  % S x S
covVm = cov(Vm');  % S x S
cCovV = bsxfun(@minus, Vm, mean(Vm,2)) * Vc' / (size(Vc, 2) - 1);  % S x S
covP = sum((U * cCovV) .* U, 2)';  % 1 x P
varP1 = sum((U * covVc) .* U, 2)';  % 1 x P
varP2 = sum((U * covVm) .* U, 2)';  % 1 x P
stdPxPy = varP1 .^ 0.5 .* varP2 .^ 0.5; % 1 x P
cMap = gather((covP ./ stdPxPy)');

% movie for predicted variance
cMovie = zeros(size(U,1),frames, 'single');
for iFrames = 1:frames
    
    frameIdx = iFrames:frames:size(Vc,2); %index for the same frame in each trial
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
    
end
fprintf('Run finished. RMSE: %f\n', median(cMovie(:).^2));


end
end