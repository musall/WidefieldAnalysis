function BpodImager_moveRawBehaviorVideo(Animal, sPath, tPath)
% Code to move raw video data off the server. 
% Finds BehaviorVideo folders in all recordings of the sourcePath and
% confirms that they contain SVD and pupil data already. If so, it moves
% the raw video to the target path and also copies frametimes and the bpod
% file.

if ~exist('sPath','var') || isempty(sPath)
    sPath = 'Y:\data\BpodImager\Animals\';
%     sPath = 'U:\smusall\BpodImager\Animals';
end
if ~exist('tPath','var') || isempty(tPath)
%     tPath = 'X:\smusall\BehaviorVideo\Animals\';
    tPath = 'Y:\TapeDrive7\BehaviorVideo\Animals';
end

if sPath(end) ~= filesep
    sPath(end + 1) = filesep;
end

if tPath(end) ~= filesep
    tPath(end + 1) = filesep;
end

%%
csPath = [sPath Animal filesep 'SpatialDisc' filesep]; %current source path
recs = dir(csPath);

for iRecs = 1 : length(recs)
    
    if strcmp(recs(iRecs).name, '.') || strcmp(recs(iRecs).name, '..')
        cFiles = [];
    else
        cFiles = dir([csPath recs(iRecs).name]); %get files/folders of current recording
    end
    
    % search for behavior video folder
    checker = false(1,3); %flag to move raw video data from current recording if SVD/pupildata is present
    for iFiles = 1 : length(cFiles)
        if strcmp(cFiles(iFiles).name, 'BehaviorVideo')
            
            cVideos = dir([csPath recs(iRecs).name filesep cFiles(iFiles).name]);
            for iVids = 1 : length(cVideos)
                if strcmpi(cVideos(iVids).name, 'motionSVD_CombinedSegments.mat')
                    checker(1) = true;
                elseif strcmpi(cVideos(iVids).name, 'SVD_CombinedSegments.mat')
                    checker(2) = true;
                elseif strcmpi(cVideos(iVids).name, 'FilteredPupil.mat')
                    checker(3) = true;
                end
            end
            
            if sum(checker) == length(checker) %if all conditions are true, move/copy some data
                
                rawVids = dir([csPath recs(iRecs).name filesep cFiles(iFiles).name filesep Animal '*mj2']);
                rawMP4s = dir([csPath recs(iRecs).name filesep cFiles(iFiles).name filesep Animal '*mp4']);
                rawVids = [rawVids; rawMP4s];
                frameTimes = dir([csPath recs(iRecs).name filesep cFiles(iFiles).name filesep Animal '*frameTimes*.mat']);
                bhvFile = dir([csPath recs(iRecs).name filesep Animal '*.mat']);
                
                if length(rawVids) > 2
                    disp('==================');
                    disp(['Found behavior SVD: ' csPath recs(iRecs).name filesep cFiles(iFiles).name])
                    disp(['Moving ' num2str(length(rawVids)) ' video files']);
                    disp(['Copying ' num2str(length(frameTimes)) ' timestamp files']);
                    disp(['Copying ' num2str(length(bhvFile)) ' bpod file(s)']);
                    tic;
                    
                    if ~exist([tPath Animal filesep recs(iRecs).name filesep], 'dir')
                        mkdir([tPath Animal filesep recs(iRecs).name filesep]);
                    end
                    
                    % show indicator positions to ensure eyeTrace ect will be correct.
                    faceVids = dir([csPath recs(iRecs).name filesep cFiles(iFiles).name filesep Animal '*1.mj2']); %cam1 should be the face camera
                    if isempty(faceVids)
                    faceVids = dir([csPath recs(iRecs).name filesep cFiles(iFiles).name filesep Animal '*1.mp4']); %cam1 should be the face camera
                    end
                    v = VideoReader([csPath recs(iRecs).name filesep cFiles(iFiles).name filesep faceVids(1).name]); %get single frame from facecam
                    pic = readFrame(v); clear v;
                    
                    load([csPath recs(iRecs).name filesep bhvFile(1).name]); %load bhv data
                    h = Behavior_checkNose(pic, SessionData); %show indicator positions
                    h.Name = recs(iRecs).name;
                    h.WindowStyle = 'docked';
                    drawnow;
                    
                    %move raw video files
                    for iVids = 1 : length(rawVids)
                        sourceFile = [csPath recs(iRecs).name filesep cFiles(iFiles).name filesep rawVids(iVids).name];
                        targetFile = [tPath Animal filesep recs(iRecs).name filesep rawVids(iVids).name];
                        if contains(sourceFile, '0001_1.mj2') || contains(sourceFile, '0001_2.mj2')
                            copyfile(sourceFile, targetFile); %dont move first two video files
                        else 
                            movefile(sourceFile, targetFile);
                        end
                    end
                    
                    %copy timestamp files
                    for iVids = 1 : length(frameTimes)
                        sourceFile = [csPath recs(iRecs).name filesep cFiles(iFiles).name filesep frameTimes(iVids).name];
                        targetFile = [tPath Animal filesep recs(iRecs).name filesep frameTimes(iVids).name];
                        copyfile(sourceFile, targetFile);
                    end
                    
                    %copy bpod file
                    for iVids = 1 : length(bhvFile)
                        sourceFile = [csPath recs(iRecs).name filesep bhvFile(iVids).name];
                        targetFile = [tPath Animal filesep recs(iRecs).name filesep bhvFile(iVids).name];
                        copyfile(sourceFile, targetFile);
                    end
                    
                    toc;
                    disp('Done');
                elseif length(rawVids) ~= 2
                    try %if no video data is found, try to copy first set of videos back to the server
                        cFile = dir([tPath Animal filesep recs(iRecs).name filesep '*0001_1.mj2']);
                        sourceFile = [tPath Animal filesep recs(iRecs).name filesep cFile.name];
                        targetFile = [csPath recs(iRecs).name filesep 'BehaviorVideo' filesep cFile.name];
                        copyfile(sourceFile, targetFile);
                        
                        cFile = dir([tPath Animal filesep recs(iRecs).name filesep '*0001_2.mj2']);
                        sourceFile = [tPath Animal filesep recs(iRecs).name filesep cFile.name];
                        targetFile = [csPath recs(iRecs).name filesep 'BehaviorVideo' filesep cFile.name];
                        copyfile(sourceFile, targetFile);
                        
                        disp('==================');
                        disp(['Found behavior SVD: ' csPath recs(iRecs).name filesep cFiles(iFiles).name])
                        disp('No movie files found. Copying two videos back to server.');
                    end
                end
            end
        end
    end
end
