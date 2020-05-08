function Widefield_batchBlockSVD(sPath, recType)

if sPath(end) ~= filesep
    sPath(end + 1) = filesep;
end

if ~exist('recType', 'var') || isempty(recType)
    recType = 'widefield';
end

nrBlocks = 49;
overlap = 10;
trigChannels = [6 7];

%%
animals = dir(sPath);

for iAnimals = 1 : length(animals)
    if ~(strcmp(animals(iAnimals).name, '.') || strcmp(animals(iAnimals).name, '..') || ~animals(iAnimals).isdir)
        
        params = dir([sPath animals(iAnimals).name filesep]);
        
        for iParams = 1 : length(params)
            if ~(strcmp(params(iParams).name, '.') || strcmp(params(iParams).name, '..') || ~params(iParams).isdir)
                
                recs = dir([sPath animals(iAnimals).name filesep params(iParams).name filesep]);
                
                for iRecs = 1 : length(recs)
                    if ~(strcmp(recs(iRecs).name, '.') || strcmp(recs(iRecs).name, '..') || ~recs(iRecs).isdir)
                    
                        cPath = [sPath animals(iAnimals).name filesep params(iParams).name filesep recs(iRecs).name filesep];
                        try
                            checker = false; %check for single or dual channel recording
                            [~,Analog] = Widefield_LoadData([cPath 'Analog_1.dat'],'Analog'); %load analog data
                            
                            if size(Analog,1) >= max(trigChannels) % assume double channel recording if both channels contain the same nr of +- 10 frames
                                checker = abs(sum(diff(Analog(trigChannels(1),:) > 1000) == 1) ...
                                          - sum(diff(Analog(trigChannels(2),:) > 1000) == 1)) <= 10; 
                            end
                                
                            if checker
                                Widefield_saveRawToBlocks(cPath, recType, nrBlocks, overlap);
                            else
                                Widefield_saveSingleChanToBlocks(cPath, recType, nrBlocks, overlap)
                            end
                        catch ME
                            disp(['SVD failed. Recording: ' cPath]);
                            disp(ME.message);
                        end
                        
                    end
                end
            end
        end
    end
end
end