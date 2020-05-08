% function BpodImager_testHemoDims
% code to test dimensionality of GLM prediction for either GFP or GCamP
% data. First test how many dimenions can be predicted in either case.

sRate = 30;
ridgeFolds = 10;

for iMod = 1:3
    
    if iMod == 1
        % GFP data
        [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
        dataOverview = dataOverview(ismember(dataOverview(:,4), 'GFP'), :);
    elseif iMod == 2
        % camK data
        [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
        dataOverview = dataOverview(ismember(dataOverview(:,4), 'gcamp'), :);
    elseif iMod == 3
        %ai93 data
        [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings;
%         cPath = 'U:\smusall\BpodImager\Animals\';
    end

    for iAnimal = 1 : size(dataOverview,1)
        %% load model data
        opts.fPath = [cPath dataOverview{iAnimal,1} filesep 'SpatialDisc' filesep dataOverview{iAnimal,3} filesep];
        disp(opts.fPath);
        load([opts.fPath 'interpVc.mat'], 'Vc');
        load([opts.fPath 'dnVc.mat'], 'U');
        load([opts.fPath 'dimBeta.mat'], 'dimBeta', 'ridgeVals');
        load([opts.fPath 'regData.mat'], 'fullR');
        load([opts.fPath 'mask.mat'], 'mask');
        randIdx = randperm(size(Vc,2)); %generate randum number index if required
        U = arrayShrink(U, mask, 'merge');
        
        %% run cross-validation
        Vm = zeros(size(Vc),'single');
        foldCnt = floor(size(Vc,2) / ridgeFolds);
        
        tic
        for iFolds = 1:ridgeFolds
            tic
            dataIdx = true(1,size(Vc,2));
            
            if ridgeFolds > 1
                
                dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
                [~, betas] = ridgeMML(Vc(:,dataIdx)', fullR(dataIdx,:), true, ridgeVals); %get beta weights for current model. Re-center and supply ridge penalty.
                Vm(:,~dataIdx) = (fullR(~dataIdx,:) * betas)'; %predict remaining data
                
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
        corrMaps{iMod}(:,:,iAnimal) = arrayShrink(cMap,mask,'split');
        avgCorr{iMod}(iAnimal) = nanmean(cMap(:).^2);
        
        [dimVar{iMod, iAnimal}, dimIdx] = sort(nanmean(Vc.^2,2),'descend');
        Vc = Vc(dimIdx,:);
        Vm = Vm(dimIdx,:);
        modVar{iMod} = nanmean(Vm.^2,2);
        
        for iDims = 1 : size(Vm,1)
            predVar{iMod}(iDims, iAnimal) = corr2(Vc(iDims,:), Vm(iDims,:))^2; %compute predicted variance for each dimension
        end        
    end
end

