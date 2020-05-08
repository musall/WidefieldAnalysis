BpodImager_checkFullModel


%% loop through data and get full model results
for iMod = 1:3
    
    if iMod == 1
        % GFP data
        [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
        dataOverview = dataOverview(ismember(dataOverview(:,4), 'GFP'), :);
    elseif iMod == 2
        % camK data
        [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
        dataOverview = dataOverview(ismember(dataOverview(:,4), 'gcamp'), :);
    elseif iMod == 3
        %ai93 data
        [dataOverview, ~, ~, ~, ~, ~, ~, ~] = delayDecRecordings;
        cPath = 'U:\smusall\BpodImager\Animals\';
    end
    
    for iAnimal = 1 : size(dataOverview,1)
        
        opts.fPath = [cPath dataOverview{iAnimal,1} filesep 'SpatialDisc' filesep dataOverview{iAnimal,3} filesep];
        disp(opts.fPath);
        
        load([opts.fPath filesep 'mask.mat']);
        load([opts.fPath filesep 'predVariance' filesep 'shCurrent' filesep 'fullCorr.mat']);
        corrMaps{iMod}(:,:,iAnimal) = arrayShrink(cMap,mask,'split');
        avgCorr{iMod}(iAnimal) = nanmean(cMap(:).^2);

        load([opts.fPath filesep 'predVariance' filesep 'shCurrent' filesep 'motorCorr.mat']);
        taskCorrMaps{iMod}(:,:,iAnimal) = arrayShrink(cMap,mask,'split');
        avgTaskCorr{iMod}(iAnimal) = nanmean(cMap(:).^2);

        load([opts.fPath filesep 'predVariance' filesep 'shOther' filesep 'motorCorr.mat']);
        motorCorrMaps{iMod}(:,:,iAnimal) = arrayShrink(cMap,mask,'split');
        avgMotorCorr{iMod}(iAnimal) = nanmean(cMap(:).^2);
        
    end
    

end

%% show some results

