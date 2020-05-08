function trials = Widefield_CheckFrameNrs(path,fileName)
%short code to identify trialnumbers from binary files in the target folder
%'path'. fileName is the name of the file before the trial number. Default
%is 'Frames_' for widefield imaging data.

if nargin < 2
    fileName = 'Frames';
end

recs = ls([path '\' fileName '*']);
temp = reshape(recs',1,numel(recs));
temp = strrep(temp,'.dat','    ');
temp = reshape(temp,size(recs'))';
if ~isempty(temp)
    trials = sort(str2num(temp(:,strfind(temp(1,:),'_')+1:end))); %identified trials
else
    trials = [];
end
