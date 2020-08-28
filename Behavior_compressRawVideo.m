
function Behavior_compressRawVideo(sPath, Animal)
% Code to compress raw video data off the server.
% Finds BehaviorVideo folders in all recordings of the sourcePath and
% confirms that they contain SVD and pupil data already. If so, it moves
% the raw video to the target path and also copies frametimes and the bpod
% file.

if ~exist('sPath','var') || isempty(sPath)
    sPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\';
end

if sPath(end) ~= filesep
    sPath(end + 1) = filesep;
end

%%
recs = dir(sPath);
for iRecs = 1 : length(recs)
    
    if strcmp(recs(iRecs).name, '.') || strcmp(recs(iRecs).name, '..')
        rawVids = [];
    else
        rawVids = dir([sPath recs(iRecs).name filesep Animal '*mj2']); %get raw videos of current recording
        fprintf('Animal: %s; Rec: %s\nCompressing %d video files\n', Animal, recs(iRecs).name, length(rawVids));
    end
    
    % search for behavior video folder
    for iVids = 1 : length(rawVids)
        try
            cFile = [sPath recs(iRecs).name filesep rawVids(iVids).name];
            rawData = squeeze(importdata(cFile));
            if max(rawData(:)) > 255 %data is not int8
                if max(rawData(:)) > 65535 %data is not int16
                    error('Unknown data type');
                else
                    rawData = mat2gray(rawData,[0 65535]);
                end
            else
                rawData = mat2gray(rawData,[0 255]);
            end
            
            rawData = padarray(rawData, 8 - [rem(size(rawData,1),8), rem(size(rawData,2),8)], 0, 'post'); %pad array to be dividable by 8
            rawData = reshape(rawData, size(rawData,1),size(rawData,2),1,[]);
            
            v = VideoWriter(strrep(cFile, '.mj2', ''), 'MPEG-4'); %save as compressed video file
            v.Quality = 100;
            open(v);
            if size(rawData,4) > 500
                for iFrames = 1 : size(rawData,4)
                    writeVideo(v,rawData(:,:,:,iFrames));
                end
            else
                writeVideo(v,rawData);
            end
            close(v);
            delete(cFile);
        end
        
    end
end