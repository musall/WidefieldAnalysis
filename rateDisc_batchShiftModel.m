function rateDisc_batchShiftModel(animal)

if ispc
    cPath = '\\grid-hs\churchland_nlsas_data\data\BpodImager\Animals\'; %data path on the server
else
    cPath = '/sonas-hs/churchland/nlsas/data/data/BpodImager/Animals/';
end

aPath = [cPath animal filesep];
recs = dir([cPath animal filesep 'SpatialDisc' filesep]);
fprintf('Basepath: %s; Found %d recordings\n', aPath, length(recs));

for iRecs = 1 : length(recs)
    
    fPath = [cPath animal filesep 'SpatialDisc' filesep recs(iRecs).name filesep];
    if exist([fPath 'interpVc.mat'], 'file') == 2
        
        fprintf('Datapath: %s\n', recs(iRecs).name);
        for iShift = 1 : 45
            
            cLine = ['qsub -l m_mem_free=1G -pe threads 8 -binding linear:8 shiftModel.sh ' ...
                animal ' ' recs(iRecs).name ' ' num2str(iShift)];
            system(cLine);
            
            cLine = ['qsub -l m_mem_free=1G -pe threads 8 -binding linear:8 shiftModel.sh ' ...
                animal ' ' recs(iRecs).name ' ' num2str(-iShift)];
            system(cLine);
            
        end
    end
end