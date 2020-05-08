function BpodImager_showPSTH(cMod)

%% select data sets
dataOverview = {
    
% 'mSM30' 'Visual' '05-Aug-2017'; ...
% 'mSM30' 'Visual' '10-Aug-2017'; ...
% 'mSM39' 'Visual' '08-Aug-2017'; ...
% 'mSM39' 'Visual' '10-Aug-2017'; ...
% 'mSM40' 'Visual' '03-Aug-2017'; ...
% 'mSM40' 'Visual' '11-Aug-2017'; ...
'mSM34' 'Visual' '31-Jul-2017'; ...
% 'mSM34' 'Visual' '21-Aug-2017'; ...

% 'mSM32' 'Audio' '17-Aug-2017'; ...
% 'mSM32' 'Audio' '18-Aug-2017'; ...
% 'mSM42' 'Audio' '29-Jul-2017'; ...
% 'mSM42' 'Audio' '31-Jul-2017'; ...
'mSM43' 'Audio' '31-Aug-2017'; ...
'mSM43' 'Audio' '12-Sep-2017'; ...
'mSM44' 'Audio' '12-Sep-2017'; ...
'mSM44' 'Audio' '31-Aug-2017'; ...

'mSM33' 'Mixed' '02-Aug-2017'; ...
'mSM33' 'Mixed' '15-Aug-2017'; ...
'mSM35' 'Mixed' '04-Aug-2017'; ...
'mSM35' 'Mixed' '07-Aug-2017'; ...
'mSM36' 'Mixed' '17-Aug-2017'; ...
'mSM36' 'Mixed' '18-Aug-2017'; ...

};

%% general variables
Paradigm = 'SpatialDisc';
cPath = 'H:\BpodImager\Animals\'; %Widefield data path
sPath = 'U:\space_managed_data\BpodImager\Animals\'; %Widefield data path on grid server
pxPerMM = 206/4;    % pixels per milimeters - default is 206/4
rotAngle = 40;      % Angle to rotate imaging data

if strcmpi(cMod,'visual')
    vars.mod = 1;
elseif strcmpi(cMod,'audio')
    vars.mod = 2;
elseif strcmpi(cMod,'mixed')
    vars.mod = 3;
end

vars.side = 1;
vars.reward = true;

%%
animals = dataOverview(:,1);
Cnt = 0;

for iAnimals = 1:length(animals)
    if strcmpi(dataOverview{iAnimals,2},cMod)
        %% check opts file, create one if not present
        Cnt = Cnt+1;
        fPath = [cPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
        disp(fPath);
        load([fPath 'Snapshot_1.mat']);
        check = false;

        if isempty(dir([fPath animals{iAnimals} '_opts.mat'])) %check if opts file exist
            check = true;
            sourcePath = [sPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep];
        else
            load([fPath animals{iAnimals}  '_opts.mat']);
            
            if ~isfield(opts,'bregma')
                check = true;
            end
        end
        
        if check
            while check
                
                opts.modality = questdlg('Select animal modality','Animal expertise','Visual','Audio','Mixed','Visual');
                opts.pxPerMM = pxPerMM;
                
                copyfile([sourcePath 'Snapshot_1.mat'],[fPath 'Snapshot_1.mat']);
                load([fPath '\Snapshot_1.mat']);
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
            save([fPath animals{iAnimals} '_opts.mat'], 'opts')
        end
        
        %% load data
        fPath = [cPath animals{iAnimals} filesep Paradigm filesep dataOverview{iAnimals,3} filesep dataOverview{iAnimals,1} '-' dataOverview{iAnimals,3}];
        load([fPath '-psthAudio.mat'])
        load([fPath '-psthVision.mat'])
        
         %load behavior data
        bhvFile = dir([fPath '*SpatialDisc*Session*.mat']);
        load([fPath bhvFile(1).name]);
        
        % get indices
        idxReward = (SessionData.Rewarded == vars.reward);
        idxSide = (SessionData.ResponseSide == vars.side);
        idxMod = (SessionData.StimType == vars.mod);
        
        %% collect data
        bhv{Cnt} = SessionData;
        allV{Cnt} = Vc;
        allU{Cnt} = U;
        aBeta{Cnt} = dimBeta;
        aMask{Cnt} = mask;
        aIdx{Cnt} = idx;
        aRecIdx{Cnt} = recIdx;
        
    end
end
        
        
snap = mapAlign(snap,opts);
 