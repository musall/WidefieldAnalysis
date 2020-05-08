function BpodImager_SVD_PSTH(Animal,Rec)
% Code to analyze SVD compressed widefield data

%% basic variables
% Animal = 'mSM39';
% Rec = '13-Jul-2017';
Paradigm = 'SpatialDisc';
cPath = ['H:\BpodImager\Animals\' Animal filesep Paradigm filesep Rec]; %Widefield data path
sPath = ['U:\space_managed_data\BpodImager\Animals\' Animal filesep Paradigm filesep Rec ]; %server data path
rotAngle = 40;      % Angle to rotate imaging data
dimCnt = 250;       % number of dimensions used for analysis
pxPerMM = 206/4;    % pixels per milimeters - default is 206/4

%% load SVD data
if isempty(dir([cPath '\Vc.mat'])) %check if file exists on hdd and pull from server otherwise
    mkdir(cPath);
    copyfile([sPath '\Vc.mat'],[cPath '\Vc.mat']);
    copyfile([sPath '\mask.mat'],[cPath '\mask.mat']);
end

disp(cPath);
load([cPath '\Vc.mat'])
load([cPath '\mask.mat'])

%% check opts file, create one if not present
check = false;
if isempty(dir([cPath '\' Animal '_opts.mat'])) %check if opts file exist
    check = true;
else
    load([cPath '\' Animal '_opts.mat']);
    
    if ~isfield(opts,'bregma')
        check = true;
    end
end

if check
    while check
        
        opts.modality = questdlg('Select animal modality','Animal expertise','Visual','Audio','Mixed','Visual');
        opts.pxPerMM = pxPerMM;
        
        copyfile([sPath '\Snapshot_1.mat'],[cPath '\Snapshot_1.mat']);
        load([cPath '\Snapshot_1.mat']);
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
    save([cPath '\' Animal '_opts.mat'], 'opts')
end

%% load behavior data and get indices for different trial types
cFile = ls([cPath filesep Animal '_' Paradigm '*.mat']);
if isempty(cFile)
    sFile = ls([sPath filesep Animal '_' Paradigm '*.mat']);
    copyfile([sPath filesep sFile],[cPath filesep sFile]);
    cFile = ls([cPath filesep Animal '_' Paradigm '*.mat']);
end
load([cPath filesep cFile]); %load behavior data

ind = trials > SessionData.nTrials; %check for imaging trials that are not in the bhv file
trials(ind) = [];

SummaryFigure(SessionData,[cPath filesep Animal '-' Rec '-BehaviorOverview.pdf'],true); %behavior summary

%% find indices
labels ={'vLeft' 'vRight' 'vSuccess' 'vError' 'aLeft' 'aRight' 'aSuccess' 'aError'};

vInd = SessionData.StimType(trials) == 1; % visual trials
aInd = SessionData.StimType(trials) == 2; % audio trials
allInd = vInd | aInd; %all unisensory trials
sucInd = SessionData.Rewarded(trials) & SessionData.Assisted(trials) & allInd; %find succesful unisensory trials
leftInd = SessionData.Assisted(trials) & ...
    ((SessionData.CorrectSide(trials) == 1 & SessionData.Rewarded(trials)) | ... %animal went left and got rewarded
    (SessionData.CorrectSide(trials) == 2 & SessionData.Punished(trials))) ; %animal went left instead of right

%% compute PSTHs
allData = zeros(size(U,1),size(U,2),size(Vc,2),length(labels),'single');

allData(:,:,:,1) = svdFrameReconstruct(U(:,:,1:dimCnt), mean(Vc(1:dimCnt,:,vInd & leftInd),3));  %average movie over all trials
allData(:,:,:,2) = svdFrameReconstruct(U(:,:,1:dimCnt), mean(Vc(1:dimCnt,:,vInd & ~leftInd),3)); %average movie over all trials
allData(:,:,:,3) = svdFrameReconstruct(U(:,:,1:dimCnt), mean(Vc(1:dimCnt,:,vInd & sucInd),3));   %average movie over all trials
allData(:,:,:,4) = svdFrameReconstruct(U(:,:,1:dimCnt), mean(Vc(1:dimCnt,:,vInd & ~sucInd),3));  %average movie over all trials

allData(:,:,:,5) = svdFrameReconstruct(U(:,:,1:dimCnt), mean(Vc(1:dimCnt,:,aInd & leftInd),3));  %average movie over all trials
allData(:,:,:,6) = svdFrameReconstruct(U(:,:,1:dimCnt), mean(Vc(1:dimCnt,:,aInd & ~leftInd),3)); %average movie over all trials
allData(:,:,:,7) = svdFrameReconstruct(U(:,:,1:dimCnt), mean(Vc(1:dimCnt,:,aInd & sucInd),3));   %average movie over all trials
allData(:,:,:,8) = svdFrameReconstruct(U(:,:,1:dimCnt), mean(Vc(1:dimCnt,:,aInd & ~sucInd),3));  %average movie over all trials

