function Widefield_checkMask(cPath, skipCheck, circleMask)
% short code to check if a 'mask.mat' file exists in a given imaging folder
% and create one if required. cPath should point directly to the folder
% that has the raw imaging data.

if ~exist('skipCheck','var') || isempty(skipCheck)
    skipCheck = false;
end

if ~exist('circleMask','var') || isempty(circleMask)
    circleMask = false;
end

if strcmpi(cPath(end),filesep)
    cPath = cPath(1:end-1);
end
opts.fPath = cPath;
opts.fName = 'Frames';
opts.preStim = 0.5;
opts.postStim = 0;
opts.trigLine = NaN;
opts.plotChans = true;
opts.maskThresh = 35;

if ~exist([cPath filesep 'mask.mat'],'file') || skipCheck
    
    % check if .dat or .mj2 files are present. Use slower splitvideo code for latter case.
    rawCheck = dir([opts.fPath filesep opts.fName '*dat']);
    vidCheck = dir([opts.fPath filesep opts.fName '*mj2']);
    analogCheck = dir([opts.fPath filesep 'Analog_*.dat']);
    
    if size(rawCheck,1) == size(analogCheck,1)
        blueData = Widefield_SplitChannels(opts,1);
        blueData = blueData(:,:,1:5);
    elseif size(vidCheck,1) == size(analogCheck,1)
        cFile = [opts.fPath filesep opts.fName '_1.mj2']; %current file to be read
        v = VideoReader(cFile);
        for x = 1:5
            blueData = readFrame(v);
        end
        clear v
    else
        error('Unequal number of imaging and analog data files. Aborted')
    end
    
    h = figure('name', cPath, 'renderer','painter');
    checker = true;
    blueData = median(single(squeeze(blueData)),3);

    if circleMask
        while checker
            mask = Widefield_CircularMask(blueData);
            imshow(mat2gray(blueData)); hold on;
            contour(mask);
            cOut = questdlg('happy?');
            if strcmp(cOut, 'Yes')
                checker = false;
            elseif strcmp(cOut, 'Cancel')
                return;
            end
        end
    else
        % make mask
        orgBlueData = blueData;
        mask = Widefield_ManualMask(blueData);
        blueData(~mask) = 0;
        
        while checker
            trace = smooth2a(double(blueData),20,20); %smoothed mean image to create mask
            mask = ~imfill(trace > prctile(trace(:),opts.maskThresh),'holes');
            
            imshow(mat2gray(orgBlueData)); caxis([0 0.5]); hold on;
            contour(mask,'r'); hold on;
            
            cOut = questdlg('happy?');
            if strcmp(cOut, 'Yes')
                checker = false;
            elseif strcmp(cOut, 'Cancel')
                return;
            else
                opts.maskThresh = str2num(char(inputdlg('Enter mask threshold (%)')));
            end
        end
    end
    save([opts.fPath filesep 'mask.mat'],'mask');
end

