function BpodImager_moveRawWidefield(Animal, sPath, tPath, nrRecs, copyOnly)
% Code to move raw widefield data off the server.
% Finds .mj2 files in all recordings of the sourcePath and
% confirms that they contain lowD data already. If so, it moves
% the raw video to the target path and also creates a copy of all other
% files in the folder.

if ~exist('sPath','var') || isempty(sPath)
    sPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';
%     sPath = 'U:\smusall\BpodImager\Animals';
end

if ~exist('tPath','var') || isempty(tPath)
%     tPath = 'X:\smusall\BpodImager\Animals\';
    tPath = '\\grid-hs\churchland_nlsas_data\TapeDrive\BpodImager\Animals\';    
end

if ~exist('copyOnly','var') || isempty(copyOnly)
    copyOnly = false; %flag to only copy raw video instead of moving it. Default is false.
end

if ~exist('nrRecs','var') || isempty(nrRecs)
    nrRecs = inf; %number of recordings to be moved. Default is all that are found.
end

if sPath(end) ~= filesep
    sPath(end + 1) = filesep;
end

if tPath(end) ~= filesep
    tPath(end + 1) = filesep;
end

vidName = 'Frames'; %name for raw widefield (mj2) files

%%
csPath = [sPath Animal filesep 'SpatialDisc' filesep]; %current source path
tsPath = [tPath Animal filesep 'SpatialDisc' filesep]; %current target path
recs = dir(csPath);
Cnt = 0; %counter for moved recordings

for iRecs = 1 : length(recs)
    
    if strcmp(recs(iRecs).name, '.') || strcmp(recs(iRecs).name, '..')
        cFiles = [];
    else
        cFiles = dir([csPath recs(iRecs).name]); %get files/folders of current recording
    end
    
    % search for lowD video file
    checker = false; %flag to move raw video data if SVD is present
    for iFiles = 1 : length(cFiles)
        if strcmpi(cFiles(iFiles).name, 'Vc.mat')
            checker = true;
        end
    end
    
    if checker && Cnt < nrRecs %if all conditions are true, move/copy some data

        rawVids = dir([csPath recs(iRecs).name filesep vidName '*mj2']); %find raw data in the folder
        nonVids = dir([csPath recs(iRecs).name filesep]); %find all files in the folder
        nonVids = nonVids(~ismember({nonVids.name}, {'.' '..'}));
        nonVids = nonVids(~ismember({nonVids.name}, {rawVids.name})); %exclude vido names from index
        
        if ~isempty(rawVids)
            Cnt = Cnt + 1;
            disp('==================');
            disp(['Found Vc.mat file: ' csPath recs(iRecs).name])
            tic;
            
            if ~exist([tsPath filesep recs(iRecs).name filesep], 'dir')
                mkdir([tsPath filesep recs(iRecs).name filesep]);
            end
            
            %move raw video files
            disp(['Moving ' num2str(length(rawVids)) ' video files']);
            for iFiles = 1 : length(rawVids)
                sourceFile = [csPath recs(iRecs).name filesep rawVids(iFiles).name];
                targetFile = [tsPath recs(iRecs).name filesep rawVids(iFiles).name];
                if copyOnly
                    copyfile(sourceFile, targetFile);
                else
                    movefile(sourceFile, targetFile);
                end
            end
            
            % copy non-raw files/folders
            disp(['Copying ' num2str(length(nonVids)) ' other files/folders']);
            for iFiles = 1 : length(nonVids)
                if ~strcmpi(nonVids(iFiles).name, 'predVariance')
                    sourceFile = [csPath recs(iRecs).name filesep nonVids(iFiles).name];
                    targetFile = [tsPath recs(iRecs).name filesep nonVids(iFiles).name];
                    copyfile(sourceFile, targetFile);
                end
            end
            
            toc;
            disp('Done');
        end
    end
end
end
