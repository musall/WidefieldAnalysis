function Widefield_SaveData(cPath,frames,frameTimes,dataType)
% short code to save imaging data from WidefieldImager code.
% cPath is the path of the file that should be opened. frames are the
% imaging data and frameTimes are the absolute timestamps for each frame.

if ~exist('dataType','var') || isempty(dataType)
    dataType = 'uint16';
end

fID = fopen(cPath,'wb');
fwrite(fID,length(frameTimes)+length(size(frames)),'double'); %write number of expected header values
fwrite(fID,frameTimes,'double'); %write absolute timestamps of each frame
fwrite(fID,size(frames),'double'); %write size of image data array
fwrite(fID,frames,dataType); %write image data
fclose(fID);