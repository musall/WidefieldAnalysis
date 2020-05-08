function BpodImagerHPC_runTestRegs
% code to run the 'Widefield_testRegs' command over all animals in the
% current dataset. Execute from bnbdev1.
% Needs updated 'delayDecRecordings' to get all
% recordings.

for iMod = 3
if iMod == 1
    % GFP data
    dataOverview = delayDecRecordings_GFP;
    dataOverview = dataOverview(ismember(dataOverview(:,4), 'GFP'), :);
elseif iMod == 2
    % camK data
    dataOverview = delayDecRecordings_GFP;
    dataOverview = dataOverview(ismember(dataOverview(:,4), 'gcamp'), :);
elseif iMod == 3
    % ai93 data
    dataOverview = delayDecRecordings;
end

% dataOverview = delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);

for iAnimals = 1 : length(animals)

    disp([animals{iAnimals} ' -- ' recs{iAnimals}]);

    cLine = ['qsub -l m_mem_free=4G -pe threads 8 -binding linear:8 modelVariance.sh ' ...
        animals{iAnimals} ' ' recs{iAnimals} ' false'];
    system(cLine);
    
    cLine = ['qsub -l m_mem_free=4G -pe threads 8 -binding linear:8 modelVariance.sh ' ...
        animals{iAnimals} ' ' recs{iAnimals} ' true'];
    system(cLine);

end
end