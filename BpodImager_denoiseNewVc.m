%% GFP data
dataOverview = delayDecRecordings_GFP;
cPath = 'Y:\data\BpodImager\Animals\';
opts.nSVD = 200;
opts.frameRate = 30;
opts.blockDims = 500; %number of dimensions from SVD per block
opts.maxLag = 5; %lag for autocorrelation
opts.autoConfidence = 0.99; %confidence for autocorrelation test
opts.autoThresh = 1.5; %threshold for autocorrelation test
opts.snrThresh = 1.6; %threshold for SNR test

for iMod = 3
    
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
        [dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings;
        cPath = 'Y:\data\BpodImager\Animals\';
    end
    
    for iAnimal = 1 : size(dataOverview,1)
        
        opts.fPath = [cPath dataOverview{iAnimal,1} filesep 'SpatialDisc' filesep dataOverview{iAnimal,3} filesep];
        disp(opts.fPath);
               
        load([opts.fPath filesep 'blueV.mat']);
        load([opts.fPath filesep 'hemoV'], 'hemoV');
        
%         [nV, rankIdx] = Widefield_denoiseVc(cat(3,blueV,hemoV), opts);
%         disp(['Rank: ' num2str(sum(rankIdx))]);
%         blueV = nV(:, :, 1 : size(blueV,3));
%         nV(:, :, 1 : size(blueV,3)) = [];
%         hemoV = nV;
%         U = U(:,:, rankIdx);
% 
%         U = U(:, :, 1:200);
%         blueV = blueV(1:200, :, :);
%         hemoV = hemoV(1:200, :, :);
        
        [Vc, regC, T, hemoVar] = Widefield_SvdHemoCorrect(U, blueV, hemoV, opts.frameRate);
        Vc = Vc(1:200, :, :);
        U = U(:, :, 1:200);
        
        save([opts.fPath filesep 'Vc.mat'],'Vc','U','blueFrametimes','Sv','totalVar','trials','bTrials','-v7.3');
        save([opts.fPath filesep 'HemoCorrection.mat'],'regC','T', 'hemoVar')
        
    end
end
