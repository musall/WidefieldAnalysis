function Widefield_checkMaskFolder(cPath, skipCheck, circleMask)

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

if ~exist('skipCheck','var') || isempty(skipCheck)
    skipCheck = false;
end

if ~exist('circleMask','var') || isempty(circleMask)
    circleMask = false;
end


%% check for subfolders
allRecs = dir(cPath);
allRecs(1:2) = [];

for iRecs = 1:length(allRecs)
    try
        Widefield_checkMask([cPath allRecs(iRecs).name], skipCheck, circleMask)
        disp(['Completed folder: ' cPath allRecs(iRecs).name]);
    catch ME
        disp(['Error in folder: ' cPath allRecs(iRecs).name]);
        disp(ME.message);
    end
end
    
%%    

