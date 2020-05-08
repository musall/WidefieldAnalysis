function Widefield_checkSize
%code to check if produced data sets have the expected size. 
% SVD code may have been compromised otherwise.

[dataOverview, ~, ~, ~, ~, ~, ~, fPath] = delayDecRecordings;
checker = false(1,size(dataOverview,1));

for iAnimals = 1 : size(dataOverview,1)
    
    cPath = [fPath dataOverview{iAnimals,1} filesep 'SpatialDisc' filesep dataOverview{iAnimals,3} filesep];
    cCheck = dir([cPath 'firstV.mat']);
    fSize = (cCheck.bytes / 1E6);
    if fSize < 300 %this file should be at least 300mb in size
        checker(iAnimals) = true;
        disp([dataOverview{iAnimals,1} ' - ' dataOverview{iAnimals,3} '; firstV is ' num2str(fSize) 'mb'])
    end
end
        
    
