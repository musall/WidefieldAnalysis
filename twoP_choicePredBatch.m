% function twoP_choicePredBatch

[dataOverview, motorLabels, ~, ~, ~, ~, ~, cPath] = twoP_delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);
visExp = ismember(dataOverview(:,2),'Visual'); %index for visual experts
opLabels = {'lGrab' 'lGrabRel' 'rGrab' 'rGrabRel' 'lLick' 'rLick'};
stepSize = 30;

% Cnt = 0; brokenRec = [];
for iAnimals = 1 : length(animals)
    tic
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep];
    disp(fPath); %current data path

    %% predict choice using re-aligned widefield
    cFile = dir([fPath filesep animals{iAnimals} '_SpatialDisc*.mat']);
    load([fPath strtrim(cFile.name)]); %load behavior data
    load([fPath 'data.mat']);
    bTrials = data.bhvTrials;
    load([fPath 'interpVc.mat']);
    load([fPath 'regData.mat']);
    
    [alignIdx, trialIdx, frames] = Widefield_getRealignment(fullR, idx, recIdx, trialIdx, recLabels, frames);
    Vc = reshape(Vc(:,alignIdx), size(Vc,1), frames, []);
    [cpRaw{iAnimals}, betaRaw{iAnimals}] = delayDec_logRegress(SessionData, Vc, [], bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize);
    
    %load motor model
    load([fPath 'motorBeta.mat'], 'motorBeta');
    load([fPath 'motorregData.mat'], 'motorR');
    motorBeta = mean(cat(3,motorBeta{:}),3);
    Vm = (motorR * motorBeta)';
    Vm = reshape(Vm(:,alignIdx), size(Vm,1), frames, []);
    [cpMotor{iAnimals}, betaMotor{iAnimals}] = delayDec_logRegress(SessionData, Vm, [], bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize);
    [cpTask{iAnimals}, betaTask{iAnimals}] = delayDec_logRegress(SessionData, Vc-Vm, [], bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize);
    
     %load spont motor model
    load([fPath 'opMotorBeta.mat'], 'opMotorBeta');
    load([fPath 'opMotorregData.mat'], 'opMotorR');
    motorBeta = mean(cat(3,opMotorBeta{:}),3);
    Vm = (opMotorR * motorBeta)';
    Vm = reshape(Vm(:,alignIdx), size(Vm,1), frames, []);
    [cpOpMotor{iAnimals}, betaOpMotor{iAnimals}] = delayDec_logRegress(SessionData, Vm, [], bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize);
     
    %load operant motor model
    load([fPath 'spontMotorBeta.mat'], 'spontMotorBeta');
    load([fPath 'spontMotorregData.mat'], 'spontMotorR');
    motorBeta = mean(cat(3,spontMotorBeta{:}),3);
    Vm = (spontMotorR * motorBeta)';
    Vm = reshape(Vm(:,alignIdx), size(Vm,1), frames, []);
    [cpSpontMotor{iAnimals}, betaSpontMotor{iAnimals}] = delayDec_logRegress(SessionData, Vm, [], bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize);
    [cpNoSpontTask{iAnimals}, betaNoSpontTask{iAnimals}] = delayDec_logRegress(SessionData, Vc-Vm, [], bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize);

    %load video only model
    load([fPath filesep 'orgregData.mat'], 'fullR', 'recLabels', 'recIdx', 'idx');
    cIdx = ismember(recIdx(~idx), find(ismember(recLabels,{'Move' 'bhvVideo'}))); %get index for task regressors
    videoR = fullR(:, cIdx)';
    videoR = reshape(videoR(:,alignIdx), size(videoR,1), frames, []);
    cpVideo{iAnimals} = delayDec_logRegress(SessionData, videoR, [], bTrials(trialIdx), strcmpi(dataOverview{iAnimals,2}, 'Visual'), stepSize);

    toc
end


%% predict choice figure
figure
cvChoice1 = cat(3,cpRaw{:});
subplot(2,2,1)
lines(1) = stdshade(squeeze(cvChoice1(:,1,:))','k',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
lines(2) = stdshade(squeeze(cvChoice1(:,2,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5);
lines(3) = stdshade(squeeze(cvChoice1(:,3,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
vline([55 81 99 114 132 162].*(1/30));
title('10x cvChoice - Raw flourescence'); legend(lines,{'All' 'Expert' 'Novice'})
xlabel('Time(s)'); ylabel('Classifier accuracy');
clear lines

cvChoice1 = cat(3,cpNoSpontTask{:});
subplot(2,2,2)
lines(1) = stdshade(squeeze(cvChoice1(:,1,:))','k',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
lines(2) = stdshade(squeeze(cvChoice1(:,2,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5);
lines(3) = stdshade(squeeze(cvChoice1(:,3,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
vline([55 81 99 114 132 162].*(1/30));
title('10x cvChoice - Spont. corrected raw'); legend(lines,{'All' 'Expert' 'Novice'})
xlabel('Time(s)'); ylabel('Classifier accuracy');
clear lines

cvChoice1 = cat(3,cpSpontMotor{:});
subplot(2,2,3)
lines(1) = stdshade(squeeze(cvChoice1(:,1,:))','k',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
lines(2) = stdshade(squeeze(cvChoice1(:,2,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5);
lines(3) = stdshade(squeeze(cvChoice1(:,3,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
vline([55 81 99 114 132 162].*(1/30));
title('10x cvChoice - Spont. movement reconstruction'); legend(lines,{'All' 'Expert' 'Novice'})
xlabel('Time(s)'); ylabel('Classifier accuracy');
clear lines

cvChoice1 = cat(3,cpOpMotor{:});
subplot(2,2,4)
lines(1) = stdshade(squeeze(cvChoice1(:,1,:))','k',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
lines(2) = stdshade(squeeze(cvChoice1(:,2,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5);
lines(3) = stdshade(squeeze(cvChoice1(:,3,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
vline([55 81 99 114 132 162].*(1/30));
title('10x cvChoice - Inst. movement reconstruction'); legend(lines,{'All' 'Expert' 'Novice'})
xlabel('Time(s)'); ylabel('Classifier accuracy');
clear lines

%% predict novice with expert trials and vice-versa
cvChoice1 = cat(3,cpRaw{:});
figure
lines(1) = stdshade(squeeze(cvChoice1(:,4,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
lines(2) = stdshade(squeeze(cvChoice1(:,5,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
vline([55 81 99 114 132 162].*(1/30))
title('Predicted choice - Cross-modal prediction'); legend(lines,{'Expert' 'Novice'})
xlabel('Time(s)'); ylabel('Classifier accuracy');

%% show results for video prediction
figure
cvChoice1 = cat(3,cpVideo{:});
lines(1) = stdshade(squeeze(cvChoice1(:,1,:))','k',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); hold on
lines(2) = stdshade(squeeze(cvChoice1(:,2,:))','r',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5);
lines(3) = stdshade(squeeze(cvChoice1(:,3,:))','b',((1:size(cvChoice1,1)) .* (stepSize/30)) - (stepSize/30)/2,0.5); axis square
ylim([0.4 1]); hline(0.5);xlim([1/3 size(cvChoice1,1).*(stepSize/30)]);
vline([55 81 99 114 132 162].*(1/30));
title('10x cvChoice - Video prediction'); legend(lines,{'All' 'Expert' 'Novice'})
xlabel('Time(s)'); ylabel('Classifier accuracy');
clear lines