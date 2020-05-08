function BpodImager_moveNonRaw(sPath, tPath, varName)
% code to move raw data from the server or back. sPath is path of the
% source folder and tPath is the path to the destionation folder.
% if varName is provided only a single file of that name will be copied.

if ~exist('varName','var')
    varName = [];
end

if ~exist('sPath','var') || isempty(sPath)
    sPath = 'W:\data\BpodImager\Animals\';
end
if ~exist('tPath','var') || isempty(tPath)
    tPath = 'Y:\simon\BpodImager\Animals\';
end

dataOverview = delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);

for iAnimals = 1 : length(animals)
    tic
    csPath = [sPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep];
    ctPath = [tPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep];
    if ~exist(ctPath,'dir')
        mkdir(ctPath);
    end
    disp(csPath); %current data path
        
    cFiles = dir(csPath);
    cFiles(1:2) = [];
    
    % move files that are not raw data or in the BehaviorVideo folder
    for iFiles = 1 : length(cFiles)
        
        cName = cFiles(iFiles).name;
        cName(length(cName) + 1 : 100) = repmat(' ', 1, 100 - length(cName)); %make sure name is long enough :)

        if ~isempty(varName)
            if strcmpi(varName, strtrim(cName))
                copyfile([csPath cFiles(iFiles).name], [ctPath cFiles(iFiles).name]); %leave a copy of some files on the server
            end
        elseif ~(strcmpi(strtrim(cName), 'BehaviorVideo') ... %not BehaviorVideo folder
                || strcmpi(cName(1:6), 'Analog') ... %not Analog file
                || strcmpi(cName(1:6), 'Frames') ... %not Frames file
                || strcmpi(cName(1:10), 'frameTimes') ... %not frameTimes file
                || strcmpi(strtrim(cName), 'handles.mat')) %not handles file
            
            if strcmpi(strtrim(cName), 'Vc.mat') || strcmpi(strtrim(cName), 'mask.mat') ...
                    || strcmpi(cName(1:8), 'Snapshot') || strcmpi(cName(1:length(animals{iAnimals})), animals{iAnimals}) ...
                    || strcmpi(cName(1:4), 'opts') || strcmpi(cName(1:5), 'opts2')
                copyfile([csPath cFiles(iFiles).name], [ctPath cFiles(iFiles).name]); %leave a copy of some files on the server
            else
                movefile([csPath cFiles(iFiles).name], [ctPath cFiles(iFiles).name]); %move file to target folder
            end
        end
    end    
end