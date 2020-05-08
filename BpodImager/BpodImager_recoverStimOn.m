function stimOn = BpodImager_recoverStimOn(cPath,sPath,trials,stimLine)
% code to read stimulus onset time from analog data files. Recording should
% be done with WidefieldImager (8/23/17).

if ~exist('stimLine','var')
    stimLine = 6;
end

for iTrials = 1:length(trials)
    if exist([cPath '\Analog_' int2str(trials(iTrials)) '.dat'],'file') ~= 2 %check if file exists on hdd and pull from server otherwise
        copyfile([sPath '\Analog_' int2str(trials(iTrials)) '.dat'],[cPath '\Analog_' int2str(trials(iTrials)) '.dat']);
    end
    
    [~,Analog] = Widefield_LoadData([cPath '\Analog_' int2str(trials(iTrials)) '.dat'],'Analog');
    cStim = find(diff(double(Analog(stimLine,:)) > 1500) == 1);
    
    if ~isempty(cStim);
        ind = find((find(diff(double(Analog(stimLine,:)) > 1500) == -1) - cStim) > 2,1); %only use triggers that are more than 2ms long
        cStim = cStim(ind) + 1;
        stimOn(iTrials) = cStim;
    end
    
    
end