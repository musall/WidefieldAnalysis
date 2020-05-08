function BpodImager_batchRun

cPath = 'Y:\data\BpodImager\Animals\';
% cPath = 'X:\smusall\BpodImager\Animals\';
% cPath = 'U:\smusall\BpodImager\Animals\';

iMod = 3;

if iMod == 1
    % GFP data
    [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
    dataOverview = dataOverview(ismember(dataOverview(:,4), 'GFP'), :);
elseif iMod == 2
    % camK data
    [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
    dataOverview = dataOverview(ismember(dataOverview(:,4), 'gcamp'), :);
    cPath = 'U:\smusall\BpodImager\Animals\'; %Widefield data path on grid server
elseif iMod == 3
    % ai93 data
    [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings;
end

% dataOverview = delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);

% Cnt = 0; brokenRec = [];
for iAnimals = 1 : length(animals)
    tic
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep];
    disp(fPath); %current data path
        
    %% run model
    disp('RunModel'); delayDec_RegressModel(cPath,animals{iAnimals},recs{iAnimals},'Widefield');
    toc
end