function BpodImager_moveForVariance(Animal, sPath, tPath, nrRecs)
% Code to move raw widefield data off the server.
% Finds .mj2 files in all recordings of the sourcePath and
% confirms that they contain lowD data already. If so, it moves
% the raw video to the target path and also creates a copy of all other
% files in the folder.

if ~exist('sPath','var') || isempty(sPath)
%     sPath = 'Y:\data\2pData\Animals\';
    sPath = 'Y:\data\BpodImager\Animals\';
%         sPath = 'X:\smusall\BpodImager\Animals\';
end
if ~exist('tPath','var') || isempty(tPath)
%     tPath = 'U:\smusall\2pData\Animals\';
    tPath = 'U:\smusall\BpodImager\Animals\';
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

%%
csPath = [sPath Animal filesep 'SpatialDisc' filesep]; %current source path
recs = dir(csPath);
Cnt = 0; %counter for moved recordings

for iRecs = 1 : length(recs)
    
    if strcmp(recs(iRecs).name, '.') || strcmp(recs(iRecs).name, '..')
        nonVids = [];
    else
        disp(['Current folder: ' csPath recs(iRecs).name])
        nonVids = dir([csPath recs(iRecs).name filesep]); %find all files in the folder
    end
    
    if ~exist([tPath Animal filesep 'SpatialDisc' filesep recs(iRecs).name filesep], 'dir') && ~isempty(nonVids)
        mkdir([tPath Animal filesep 'SpatialDisc' filesep recs(iRecs).name filesep]);
    end
    
    %copy files
    for iFiles = 1 : length(nonVids)
        sourceFile = [csPath recs(iRecs).name filesep nonVids(iFiles).name];
        targetFile = [tPath Animal filesep 'SpatialDisc' filesep recs(iRecs).name filesep nonVids(iFiles).name];
        
        [~, b] = fileparts(sourceFile);
        
        if strcmpi(b, 'interpVc') || strcmpi(b, 'dimBeta') ||  strcmpi(b, 'orgdimBeta') ...
                || strcmpi(b, 'regData') || strcmpi(b, 'orgregData') || strcmpi(b, 'mask') || strcmpi(b, 'data') ...
                || strcmpi(b, 'Vc') || strcmpi(b, 'opts2') || strcmpi(b, 'Snapshot_1') || strcmpi(b, 'HemoCorrection')
            
            copyfile(sourceFile, targetFile);
            
        end
    end
    disp('Done');
end

end
