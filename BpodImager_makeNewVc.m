%% GFP data
dataOverview = delayDecRecordings_GFP;
% cPath = 'Y:\data\BpodImager\Animals\';
cPath = 'U:\smusall\BpodImager\Animals\';
opts.nSVD = 200;
opts.frameRate = 30;

for iAnimal = 1 : size(dataOverview,1)

opts.fPath = [cPath dataOverview{iAnimal,1} filesep 'SpatialDisc' filesep dataOverview{iAnimal,3} filesep];
disp(opts.fPath);

load([opts.fPath filesep 'blueV']);
load([opts.fPath filesep 'hemoV'], 'hemoV');
blueV = blueV(1:opts.nSVD, :, :);
hemoV = hemoV(1:opts.nSVD, :, :);
U = U(:,:, 1:opts.nSVD);
[Vc, regC, T, hemoVar] = Widefield_SvdHemoCorrect(U, blueV, hemoV, opts.frameRate);
save([opts.fPath filesep 'Vc.mat'],'Vc','U','blueFrametimes','Sv','totalVar','trials','bTrials','-v7.3');
save([opts.fPath filesep 'hemoCorrection.mat'],'regC','T', 'hemoVar')

end

%% WF data
dataOverview = delayDecRecordings;
cPath = 'X:\smusall\BpodImager\Animals\';
opts.nSVD = 200;
opts.frameRate = 30;

for iAnimal = 1 : size(dataOverview,1)
    
opts.fPath = [cPath dataOverview{iAnimal,1} filesep 'SpatialDisc' filesep dataOverview{iAnimal,3} filesep];

load([opts.fPath filesep 'blueV']);
load([opts.fPath filesep 'hemoV'], 'hemoV');
blueV = blueV(1:opts.nSVD, :, :);
hemoV = hemoV(1:opts.nSVD, :, :);
U = U(:,:, 1:opts.nSVD);
[Vc, regC, T, hemoVar] = Widefield_SvdHemoCorrect(U, blueV, hemoV, opts.frameRate);
save([opts.fPath filesep 'Vc.mat'],'Vc','U','blueFrametimes','Sv','totalVar','trials','bTrials','-v7.3');
save([opts.fPath filesep 'hemoCorrection.mat'],'regC','T', 'hemoVar')

end


