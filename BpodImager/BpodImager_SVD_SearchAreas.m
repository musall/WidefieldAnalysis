function BpodImager_SVD_SearchAreas(Animal,Rec,reload)
% Code to analyze SVD compressed widefield data. This code is to create
% PSTHs from visual and auditory data and compute dPrime maps for
% differences between modalities, choices and success. 

if ~exist('reload', 'var') || isempty(reload)
    reload = false; %flag to reload all data from server
end

%% basic variables
Paradigm = 'SpatialDisc';
% cPath = ['U:\space_managed_data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec]; %Widefield data path
cPath = ['H:\BpodImager\Animals\' Animal filesep Paradigm filesep Rec filesep]; %Widefield data path
sPath = ['U:\space_managed_data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec filesep]; %server data path
smth = 2;           %2-D smooth for plots/movies
rotAngle = 40;      % Angle to initially rotate snapshot for alignment
dimCnt = 200;       % number of dimensions used for analysis
tBin = 1;           % timebins for analysis. Used to average frames together for speed.
pxPerMM = 206/4;    % pixels per milimeters - default is 206/4

%% load SVD data
if isempty(dir([cPath 'Vc.mat'])) || reload %check if file exists on hdd and pull from server otherwise
    mkdir(cPath);
    copyfile([sPath 'Vc.mat'],[cPath 'Vc.mat']);
    copyfile([sPath 'interpVc.mat'],[cPath 'interpVc.mat']);
    copyfile([sPath 'mask.mat'],[cPath 'mask.mat']);
    copyfile([sPath 'Snapshot_1.mat'],[cPath 'Snapshot_1.mat']);
end

disp(cPath);
load([cPath 'Vc.mat'],'U','bTrials','trials')
load([cPath 'interpVc.mat'],'Vc','frames')
load([cPath 'mask.mat'])
load([cPath 'Snapshot_1.mat']);
Vc = reshape(Vc,size(Vc,1),frames,[]);

%% check for model output file
if exist([cPath 'regData.mat'],'file') ~= 2 || exist([cPath 'dimBeta.mat'],'file') ~= 2 || reload %check if file exists on hdd and pull from server otherwise
    copyfile([sPath 'regData.mat'],[cPath 'regData.mat']);
    copyfile([sPath 'dimBeta.mat'],[cPath 'dimBeta.mat']);
end

if exist([sPath Animal '_opts.mat'],'file') == 2 %if opts file exist on the server, copy to local path
    copyfile([sPath Animal '_opts.mat'],[cPath Animal '_opts.mat']);
end

load([cPath 'regData.mat'],'fullR','trialIdx','recIdx','idx','recLabels'); %load model data
load([cPath 'dimBeta.mat'],'dimBeta');
trialIdx = unique(ceil(find(trialIdx)/frames)); %convert trialIdx
fprintf('Rejecting %d/%d trials for not being used in the linear model. %d trials remaining.\n',length(trialIdx),length(bTrials),length(bTrials) - length(trialIdx));
bTrials(trialIdx) = []; %don't use trials that were rejected by the model

%% check opts file, create one if not present
check = false;
if  exist([cPath Animal '_opts.mat'],'file') ~= 2 || reload %check if opts file exist
    check = true;
else
    load([cPath Animal '_opts.mat']);
    if ~isfield(opts,'bregma')
        check = true;
    end
end

if check
    while check
        
        opts.modality = questdlg('Select animal modality','Animal expertise','Visual','Audio','Mixed','Visual');
        opts.pxPerMM = pxPerMM;
        
        copyfile([sPath 'Snapshot_1.mat'],[cPath 'Snapshot_1.mat']);
        load([cPath 'Snapshot_1.mat']);
        snap = double(imrotate(snap,rotAngle,'crop'));
        
        figure(90);
        imagesc(snap);axis image; colormap gray; hold on
        
        [bregma(1), bregma(2)] = Widefield_SelectCoordinates(snap,pxPerMM,'bregma'); %set bregma coordinate
        plot(bregma(1), bregma(2), 'or', 'MarkerFaceColor', 'r')
        snap(round(bregma(2)), round(bregma(1))) = inf;
        
        [lambda(1), lambda(2)] = Widefield_SelectCoordinates(snap,pxPerMM,'lambda'); %set lambda coordinate
        plot(lambda(1), lambda(2), 'ob', 'MarkerFaceColor', 'b')
        snap(round(lambda(2)), round(lambda(1))) = inf;
        
        [frontal(1), frontal(2)] = Widefield_SelectCoordinates(snap,pxPerMM,'frontal'); %set that other coordinate
        plot(frontal(1), frontal(2), 'ow', 'MarkerFaceColor', 'w')
        snap(round(frontal(2)), round(frontal(1))) = inf;
        
        degAngle(1) = (atan2(bregma(2) - lambda(2), bregma(1) - lambda(1))*180/pi) * -1;    %angle for lambda to bregma line, relative to x-axis
        degAngle(2) = (atan2(frontal(2) - bregma(2), frontal(1) - bregma(1))*180/pi) * -1;  %angle for bregma to frontal line, relative to x-axis
        imAngle = 90 - mean(degAngle);
        
        snap = imrotate(snap,imAngle,'crop');
        opts.rotAngle = rotAngle + imAngle; %rotation angle of raw data to straighten image
        
        [b, a] = ind2sub(size(snap),find(isinf(snap))); %find landmark coordinates from rotated image
        [b,yy] = sort(b);
        a = a(yy);
        
        opts.frontal = [a(1) b(1)];
        opts.bregma = [a(2) b(2)];
        opts.lambda = [a(3) b(3)];
        
        figure(90); hold off;
        imagesc(snap);axis image; colormap gray; hold on
        grid(gca,'on');grid minor;set(gca,'GridColor','w');
        plot(opts.bregma(1), opts.bregma(2), 'ob', 'MarkerFaceColor', 'b');
        plot(opts.lambda(1), opts.lambda(2), 'or', 'MarkerFaceColor', 'r');
        plot(opts.frontal(1), opts.frontal(2), 'ow', 'MarkerFaceColor', 'w')
        title('Aligned image');
        
        happy = questdlg('Happy?','Happy?','Yes','No','No');
        
        if strcmpi(happy,'Yes')
            check = false;
        end
    end
    save([cPath Animal '_opts.mat'], 'opts');
    save([sPath Animal '_opts.mat'], 'opts');
end

U = Widefield_mapAlign(U,opts); %rotate for straight midline and set bregma to center
snap = Widefield_mapAlign(snap,opts); %rotate for straight midline and set bregma to center
newMask = isnan(U(:,:,1)); %get new mask for aligned U

%% load behavior data and get indices for different trial types
cFile = ls([cPath Animal '_' Paradigm '*.mat']);
if isempty(cFile)
    sFile = ls([sPath filesep Animal '_' Paradigm '*.mat']);
    copyfile([sPath filesep sFile],[cPath sFile]);
    cFile = ls([cPath Animal '_' Paradigm '*.mat']);
end
load([cPath cFile]); %load behavior data

ind = bTrials > SessionData.nTrials; %check for imaging trials that are not in the bhv file
bTrials(ind) = [];

SummaryFigure(SessionData,[cPath Animal '-' Rec '-BehaviorOverview.pdf'],true); %behavior summary

%% compute delta F/F movies
trace = svdFrameReconstruct(U(:,:,1:dimCnt), mean(Vc(1:dimCnt,:,:),3));  %average movie over all trials
Widefield_SaveToAvi(trace,[Animal '-' Rec '-AllAvg-fps10'],10,'colormap_blueblackred',[-0.01 0.01],cPath,smth,2);

cMod = {'lVision' 'lAudio' 'lAudioVisual';'rVision' 'rAudio' 'rAudioVisual'};
for iMod = 1:2 % StimType = 1 (Vision), StimType = 2 (Audio), StimType = 3 (Audiovisual),
    for iSide = 1:2
        
        detectInd = SessionData.StimType(bTrials) == iMod & ~SessionData.DidNotChoose(bTrials) ...
            & ~SessionData.DidNotLever(bTrials) & SessionData.CorrectSide(bTrials) == iSide;
        trace = svdFrameReconstruct(U(:,:,1:dimCnt),median(Vc(:,:,detectInd),3)); %get average data
        trace = bsxfun(@minus,trace,  mean(trace(:,:,1:20),3)); %correct for baseline
        Widefield_SaveToAvi(trace,[Animal '-' Rec '-' cMod{iSide,iMod} '-fps10'],10,'colormap_blueblackred',[-0.01 0.01],cPath,smth,2);

    end
end

%% select trials for AUC
selTrials = Behavior_AnalysisTrialselection(SessionData,bTrials,1,2,true); %select for even number of singe modality trials and choices for each side
selInd = false(1,length(bTrials));
selInd(selTrials) = true;

selErrorTrials = Behavior_AnalysisTrialselection(SessionData,bTrials,1,2,false); %select for even number of singe modality trials and choices for each side
selErrorInd = false(1,length(bTrials));
selErrorInd(selErrorTrials) = true;

%% get trial indices
vInd = SessionData.StimType(bTrials) == 1; % visual trials
aInd = SessionData.StimType(bTrials) == 2; % audio trials
allInd = vInd | aInd; %all unisensory trials
sucInd = SessionData.Rewarded(bTrials) & SessionData.Assisted(bTrials) & allInd; %find succesful unisensory trials
leftInd = SessionData.Assisted(bTrials) & ...
    ((SessionData.CorrectSide(bTrials) == 1 & SessionData.Rewarded(bTrials)) | ... %animal went left and got rewarded
    (SessionData.CorrectSide(bTrials) == 2 & SessionData.Punished(bTrials))) ; %animal went left instead of right

save([cPath Animal '-' Rec '-trialIdx.mat'],'vInd','aInd','allInd','sucInd','selInd','selErrorInd','leftInd');

%% compute averages and save   
%vision
visAll = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,vInd),3))),newMask); %all vision trials
visCorrect = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,sucInd & vInd),3))),newMask); %all correct vision trials
visCorrectLeft = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,sucInd & vInd & leftInd),3))),newMask); %all left correct vision trials
visCorrectRight = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,sucInd & vInd & ~leftInd),3))),newMask); %all right correct vision trials
visError = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,~sucInd & vInd),3))),newMask); %all error vision trials
visErrorLeft = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,~sucInd & vInd & leftInd),3))),newMask); %all left error vision trials
visErrorRight = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,~sucInd & vInd & ~leftInd),3))),newMask); %all right error vision trials
visSelect = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,selInd & vInd),3))),newMask); %selected visual trias
save([cPath Animal '-' Rec '-psthVision.mat'],'visAll','visCorrect','visCorrectLeft','visCorrectRight','visError','visErrorLeft','visErrorRight','visSelect','newMask','snap','-v7.3');

%audio 
audAll = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,aInd),3))),newMask); %all audio trials
audCorrect = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,sucInd & aInd),3))),newMask); %all correct audio trials
audCorrectLeft = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,sucInd & aInd & leftInd),3))),newMask); %all left correct audio trials
audCorrectRight = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,sucInd & aInd & ~leftInd),3))),newMask); %all right correct audio trials
audError = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,~sucInd & aInd),3))),newMask); %all error audio trials
audErrorLeft = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,~sucInd & aInd & leftInd),3))),newMask); %all left error audio trials
audErrorRight = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,~sucInd & aInd & ~leftInd),3))),newMask); %all right error audio trials
audSelect = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vc(1:dimCnt,:,selInd & aInd),3))),newMask); %selected audio trias
save([cPath Animal '-' Rec '-psthAudio.mat'],'audAll','audCorrect','audCorrectLeft','audCorrectRight','audError','audErrorLeft','audErrorRight','audSelect','newMask','snap','-v7.3');

%% do reconstruction for all regressors
% recLabels{end+1} = 'all';
% for iRegs = 1 : length(recLabels)
%     
%     if strcmpi(recLabels{iRegs},'all')
%         cInd = true(1,size(dimBeta,1));
%     else
%         cInd = ismember(recIdx(~idx), find(ismember(recLabels,recLabels(iRegs)))); %find all motor regressors
%     end
%     
%     cData = fullR(:, cInd); %motor regressors
%     cBeta = dimBeta(cInd,:); %motor regressors
%     Vm = (cData * cBeta)'; %model Vc data
%     Vm = reshape(Vm,size(Vm,1),size(Vc,2),size(Vc,3));
%     
%     % compute averages and save
%     %vision
%     visAll = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,vInd),3))),newMask); %all vision trials
%     visCorrect = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,sucInd & vInd),3))),newMask); %all correct vision trials
%     visCorrectLeft = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,sucInd & vInd & leftInd),3))),newMask); %all left correct vision trials
%     visCorrectRight = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,sucInd & vInd & ~leftInd),3))),newMask); %all right correct vision trials
%     visError = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,~sucInd & vInd),3))),newMask); %all error vision trials
%     visErrorLeft = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,~sucInd & vInd & leftInd),3))),newMask); %all left error vision trials
%     visErrorRight = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,~sucInd & vInd & ~leftInd),3))),newMask); %all right error vision trials
%     visSelect = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,selInd & vInd),3))),newMask); %selected visual trias
%     save([cPath Animal '-' Rec '-' recLabels{iRegs} 'Rebuild-psthVision.mat'],'visAll','visCorrect','visCorrectLeft','visCorrectRight','visError','visErrorLeft','visErrorRight','visSelect','newMask','snap','-v7.3');
%     
%     %audio
%     audAll = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,aInd),3))),newMask); %all audio trials
%     audCorrect = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,sucInd & aInd),3))),newMask); %all correct audio trials
%     audCorrectLeft = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,sucInd & aInd & leftInd),3))),newMask); %all left correct audio trials
%     audCorrectRight = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,sucInd & aInd & ~leftInd),3))),newMask); %all right correct audio trials
%     audError = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,~sucInd & aInd),3))),newMask); %all error audio trials
%     audErrorLeft = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,~sucInd & aInd & leftInd),3))),newMask); %all left error audio trials
%     audErrorRight = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,~sucInd & aInd & ~leftInd),3))),newMask); %all right error audio trials
%     audSelect = arrayShrink(svdFrameReconstruct(U(:,:,1:dimCnt),squeeze(mean(Vm(1:dimCnt,:,selInd & aInd),3))),newMask); %selected audio trias
%     save([cPath Animal '-' Rec '-' recLabels{iRegs} 'Rebuild-psthAudio.mat'],'audAll','audCorrect','audCorrectLeft','audCorrectRight','audError','audErrorLeft','audErrorRight','audSelect','newMask','snap','-v7.3');
% 
% end

%% AUC analysis - Modality differences
allPrime = NaN([sum(~newMask(:)),frames/tBin]); %dPrime between all visual / audio trials
corrPrime = NaN([sum(~newMask(:)),frames/tBin]); %dPrime between all correct visual / audio trials
selPrime = NaN([sum(~newMask(:)),frames/tBin]); %dPrime between all selected visual / audio trials (correct, same nr of left/right choices in each modality)
selErrorPrime = NaN([sum(~newMask(:)),frames/tBin]); %dPrime between all selected visual / audio trials (errors, same nr of left/right choices in each modality)

visChoicePrime = NaN([sum(~newMask(:)),frames/tBin]); %dPrime between left and right in all visual trials
visSuccPrime = NaN([sum(~newMask(:)),frames/tBin]); %dPrime between correct and incorrect in all visual trials
audChoicePrime = NaN([sum(~newMask(:)),frames/tBin]); %dPrime between left and right in all audio trials
audSuccPrime = NaN([sum(~newMask(:)),frames/tBin]); %dPrime between correct and incorrect in all audio trials

U = arrayShrink(U,newMask);

% compute dPrime
Cnt = 0;
for x = 1 : tBin : frames
    
    Cnt = Cnt+1;
    temp = U(:,1:dimCnt)*squeeze(mean(Vc(1:dimCnt,x : x + tBin-1, :),2));
    
    allPrime(:,Cnt) = (mean(temp(:,allInd & vInd),2) - mean(temp(:,allInd & aInd),2)) ./ std(temp(:,allInd), [],2);
    corrPrime(:,Cnt) = (mean(temp(:,sucInd & vInd),2) - mean(temp(:,sucInd & aInd),2)) ./ std(temp(:,sucInd), [],2);
    selPrime(:,Cnt) = (mean(temp(:,selInd & vInd),2) - mean(temp(:,selInd & aInd),2)) ./ std(temp(:,selInd), [],2);

    visChoicePrime(:,Cnt) = (mean(temp(:,vInd & leftInd),2) - mean(temp(:,vInd & ~leftInd),2)) ./ std(temp(:,vInd), [],2);
    visSuccPrime(:,Cnt) = (mean(temp(:,vInd & sucInd),2) - mean(temp(:,vInd & ~sucInd),2)) ./ std(temp(:,vInd), [],2);

    audChoicePrime(:,Cnt) = (mean(temp(:,aInd & leftInd),2) - mean(temp(:,aInd & ~leftInd),2)) ./ std(temp(:,aInd), [],2);
    audSuccPrime(:,Cnt) = (mean(temp(:,aInd & sucInd),2) - mean(temp(:,aInd & ~sucInd),2)) ./ std(temp(:,aInd), [],2);
    
    if sum(selErrorInd) > 0
        selErrorPrime(:,Cnt) = (mean(temp(:,selErrorInd & vInd),2) - mean(temp(:,selErrorInd & aInd),2)) ./ std(temp(:,selInd), [],2);
    end
    
    if rem(x,tBin*10) == 0
        fprintf(1, 'dPrime: Current bin is %d out of %d\n', x, frames);
    end
end
save([cPath Animal '-' Rec '-AVmodPrime.mat'],'allPrime','corrPrime','selPrime','selErrorPrime','visChoicePrime','visSuccPrime','audChoicePrime','audSuccPrime','newMask','snap','-v7.3');


