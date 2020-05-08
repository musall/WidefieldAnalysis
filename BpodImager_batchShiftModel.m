function BpodImager_batchShiftModel

if ispc
    cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
else
    cPath = '/sonas-hs/churchland/hpc/home/smusall/BpodImager/Animals/';
%     cPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/';
end

dataOverview = delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);
    
for iRecs = 1 : length(recs)

    
    fPath = [cPath animals{iRecs} filesep 'SpatialDisc' filesep recs{iRecs} filesep];
    if exist([fPath 'interpVc.mat'], 'file') == 2
        
        fprintf('Datapath: %s\n', recs{iRecs});
        
        cLine = ['qsub -l m_mem_free=1G -pe threads 8 -binding linear:8 shiftModel.sh ' ...
            animals{iRecs} ' ' recs{iRecs} ' ' num2str(0)];
        system(cLine);
        
        for iShift = 1 : 45
            
            cLine = ['qsub -l m_mem_free=1G -pe threads 8 -binding linear:8 shiftModel.sh ' ...
                animals{iRecs} ' ' recs{iRecs} ' ' num2str(iShift)];
            system(cLine);
            
            cLine = ['qsub -l m_mem_free=1G -pe threads 8 -binding linear:8 shiftModel.sh ' ...
                animals{iRecs} ' ' recs{iRecs} ' ' num2str(-iShift)];
            system(cLine);
            
        end
    end
end