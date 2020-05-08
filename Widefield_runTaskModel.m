function Widefield_runTaskModel(Vm, cPath)

load([cPath 'regData.mat'],'fullR','recIdx','recLabels','trialIdx','idx')    
load([dataPath{iAnimals} 'interpVmotor.mat'],'Vm')

[~, motorLabels] = delayDecRecordings;
cIdx = ismember(recIdx(~idx), find(~ismember(recLabels,motorLabels))); %get index for task regressors
taskLabels = recLabels(sort(find(~ismember(recLabels,motorLabels)))); %make sure taskLabels is in the right order

%create new regressor index that matches task labels
taskIdx = recIdx(~idx);
taskIdx = taskIdx(cIdx);
temp = unique(taskIdx);
for x = 1 : length(temp)
    taskIdx(taskIdx == temp(x)) = x;
end
taskR = fullR(:,cIdx);

[taskRidge, taskBeta{iFolds}] = ridgeMML(Vm', taskR, true); %get beta weights for motor-only model.
fprintf('Mean ridge penalty for task-only, zero-mean model: %f\n', mean(taskRidge));
save([cPath 'taskdimBeta.mat'], 'taskBeta', 'taskRidge');
save([cPath 'taskregData.mat'], 'taskR','trialIdx', 'taskIdx', 'taskLabels','gaussShift','-v7.3');
    