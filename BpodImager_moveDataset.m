sPath = 'U:\space_managed_data\BpodImager\Animals\'; %sourcepath
tPath = 'F:\WidefieldData\Animals\'; %target path
dataOverview = delayDecRecordings;


for iAnimals = 1 : size(dataOverview,1)

    fPath = [sPath dataOverview{iAnimals,1} filesep 'SpatialDisc' filesep dataOverview{iAnimals,3} filesep]; %current source datapath
    tfPath = [tPath dataOverview{iAnimals,1} filesep 'SpatialDisc' filesep dataOverview{iAnimals,3} filesep]; %current target datapath
    
    if ~exist(tfPath,'dir')
        mkdir(tfPath);
    end
    
    bhvFile = strsplit(fPath,filesep);
    bhvFile = dir([fPath bhvFile{end-3} '*' bhvFile{end-2} '*.mat']);
    
    %copy some data
    copyfile([fPath 'Snapshot_1.mat'],[tfPath 'Snapshot_1.mat']);
    copyfile([fPath 'mask.mat'],[tfPath 'mask.mat']);
    copyfile([fPath 'opts2.mat'],[tfPath 'opts2.mat']);
    copyfile([fPath 'Vc.mat'],[tfPath 'Vc.mat']);
    copyfile([fPath bhvFile.name],[tfPath bhvFile.name]);
  
    disp('Finished copying data')
    disp(['Source: ' fPath]);
    disp(['Target: ' fPath]);
    disp('=====================');
    
end
