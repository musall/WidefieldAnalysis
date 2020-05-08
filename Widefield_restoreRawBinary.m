function Widefield_restoreRawBinary(cPath)
% short code to re-create binary data file from mj2 movies. This is the
% counterpiece of the Widefield_ArchiveRawData code.

fileName = 'Frames';
vidType = 'mp4';
rawVids = dir([cPath filesep fileName '*' vidType]); %get raw videos of current recording
fprintf('Rec: %s\nRebuilding %d video files\n', cPath, length(rawVids));
        
for x = 1 : length(rawVids)
    
    clear rawData frameTimes
    cFile = [cPath filesep rawVids(x).name];
    rawData = squeeze(importdata(cFile));
    
    cFile = [cPath filesep strrep(rawVids(x).name,fileName,'frameTimes')];
    cFile = strrep(cFile, vidType, 'mat');
    load(cFile);
    
    cFile = [cPath filesep rawVids(x).name];
    cFile = strrep(cFile, vidType, 'dat');

    sID = fopen(cFile, 'Wb'); %open binary stimulus file
    fwrite(sID,length(frameTimes)+length(size(rawData)),'double'); %write number of expected header values
    fwrite(sID,frameTimes,'double'); %write absolute timestamps of each frame
    fwrite(sID,size(rawData),'double'); %write size of image data array
    fwrite(sID,rawData,'uint16'); %write image data
    fclose(sID);
    disp(['done - ' cFile]);
    
end