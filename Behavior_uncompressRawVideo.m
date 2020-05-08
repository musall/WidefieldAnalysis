function Behavior_uncompressRawVideo(Animal, sPath)
% Code to compress mp4 videos back to mj2 to run dimensionality reduction
% on the HPC. Finds BehaviorVideo folders in all recordings of the sourcePath and
% confirms that they contain SVD and pupil data already. If so, it moves
% the raw video to the target path and also copies frametimes and the bpod
% file.

if ~exist('sPath','var') || isempty(sPath)
    sPath = '\\grid-hs\churchland_nlsas_data\BehaviorVideo\';
end

if sPath(end) ~= filesep
    sPath(end + 1) = filesep;
end
sPath = [sPath Animal filesep 'SpatialDisc' filesep 'Session Data' filesep];

%%
recs = dir(sPath);
for iRecs = 1 : length(recs)
    
    if strcmp(recs(iRecs).name, '.') || strcmp(recs(iRecs).name, '..')
        rawVids = [];
        svdCheck = [];
    else
        svdCheck = dir([sPath recs(iRecs).name filesep '*SVD_CombinedSegments.mat']);
        rawVids = dir([sPath recs(iRecs).name filesep Animal '*mp4']); %get raw videos of current recording
        
        if (length(svdCheck) == 2)
            disp(['Found SVDcomplete in experiment: ' recs(iRecs).name ' - skipped']);
            rawVids = [];
        else
            fprintf('Animal: %s; Rec: %s\nUncompressing %d video files\n', Animal, recs(iRecs).name, length(rawVids));
        end
    end
    
    % search for behavior video folder
    for iVids = 1 : length(rawVids)
            
        cFile = [sPath recs(iRecs).name filesep rawVids(iVids).name];
        rawData = squeeze(importdata(cFile));
        rawData = rawData(:,:,1,:);
        
        if iVids == 1
            xIdx = sum(rawData(:,:,1,1) == 0,1) > (size(rawData,1) - 20); %remove 0s if needed
            yIdx = sum(rawData(:,:,1,1) == 0,2) > (size(rawData,2) - 20); %remove 0s if needed
        end
        rawData = rawData(~yIdx,~xIdx,:,:);
        
        v = VideoWriter(strrep(cFile, '.mp4', ''), 'Archival'); %save as mj2 file
        v.FrameRate = 30;
        open(v);writeVideo(v,rawData);close(v);
        delete(cFile);
        
    end
end