% function BpodImager_crossValFullModel

dataOverview = delayDecRecordings;
animals = dataOverview(:,1)';
Paradigm = 'SpatialDisc';
cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
ridgeFolds = 10;    %number of cross-validations

%%
for iAnimals = 1 : size(dataOverview,1)

    fPath = [cPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
    load([fPath 'interpVc.mat'], 'Vc', 'frames'); %load real data
    load([fPath 'regData.mat'], 'fullR'); %load model
    load([fPath 'dimBeta.mat'], 'dimBeta'); %load weights
    
    
    % create cross-validated prediction
    Vm = zeros(size(Vc),'single'); %pre-allocate motor-reconstructed V
    randIdx = randperm(size(Vc,2)); %generate randum number index
    foldCnt = floor(size(Vc,2) / ridgeFolds);
    cBeta = cell(1,ridgeFolds);
    for iFolds = 1:ridgeFolds
        dataIdx = true(1,size(Vc,2));
        
        dataIdx(randIdx(((iFolds - 1)*foldCnt) + (1:foldCnt))) = false; %index for training data
        if iFolds == 1
            [cRidge, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', fullR(dataIdx,:), true); %get beta weights and ridge penalty for task only model
        else
            [~, cBeta{iFolds}] = ridgeMML(Vc(:,dataIdx)', fullR(dataIdx,:), true, cRidge); %get beta weights for task only model. ridge value should be the same as in the first run.
        end
        Vm(:,~dataIdx) = (fullR(~dataIdx,:) * cBeta{iFolds})'; %predict remaining data
        
        if rem(iFolds,ridgeFolds/5) == 0
            fprintf(1, 'Current fold is %d out of %d\n', iFolds, ridgeFolds);
        end
    end
    
    save([fPath 'interpVfull.mat'], 'Vm', 'cBeta', 'frames'); %save cross-validation result
    fprintf('Done. %d/%d animals\n',iAnimals,size(dataOverview,1));
end