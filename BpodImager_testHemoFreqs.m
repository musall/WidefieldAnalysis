function BpodImager_testHemoFreqs
% code to test dimensionality of GLM prediction for either GFP or GCamP
% data. First test how many dimenions can be predicted in either case.

cPath = 'X:\smusall\BpodImager\Animals\';

dataOverview = delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);

% Cnt = 0; brokenRec = [];
for iAnimals = 1 : length(animals)
    tic
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep];
    disp(fPath); %current data path
        
    disp('RunModel'); delayDec_RegressModel(cPath,animals{iAnimals},recs{iAnimals},'Widefield');
        
%         disp('BodyVars'); Behavior_computeBodyVars(cPath, animals{iAnimals}, recs{iAnimals}, 2)
%         disp('PupileVars'); Behavior_computePupilVars([fPath 'BehaviorVideo'],false)
%         Behavior_computeFaceVars(cPath, animals{iAnimals}, recs{iAnimals}, 1, false)
%         disp('RunModel'); BpodImager_delayRegressModel(cPath,animals{iAnimals},recs{iAnimals});
%         Cnt = Cnt+1;
%         brokenRec(Cnt) = iAnimals;
%     end
    toc
end



% % code to test at which frequencies the GLM can predict either GFP or GCaMP
% % data by filtering Vc at different frequencies before computing VM and
% % testing pred.Variance with full model.