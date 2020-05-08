function Widefield_mergeSplitV(cPath, rec1, rec2, rec3)

cPath = 'Y:\data\BpodImager\Animals\mSM64\SpatialDisc\';
rec1 = '10-Aug-2018';
rec2 = '10-Aug-2018_1';
rec3 = '10-Aug-2018_2';

opts.nSVD = 500;
fPath1 = [cPath rec1];
fPath2 = [cPath rec2];
newPath = [cPath rec3];

data1 = load([fPath1 filesep 'blueV.mat']);
data2 = load([fPath2 filesep 'blueV.mat']);

blueFrametimes = cat(2,data1.blueFrametimes, data2.blueFrametimes);
bTrials = [data1.bTrials data2.bTrials];
trials = [data1.trials data2.trials+data1.trials(end)];
mask1 = isnan(data1.U(:,:,1));
mask2 = isnan(data2.U(:,:,1));
mask = mask1 | mask2;


% kill later
data1.blueV = data1.blueV(:, :, 1:2);
data2.blueV = data2.blueV(:, :, 1:2);

dSize = size(data1.blueV);
data1.blueV = reshape(data1.blueV, size(data1.blueV,1), []);
data1.U = arrayShrink(data1.U, mask);
trialIdx1 = isnan(data1.blueV(1, :));
data1.blueV(:, trialIdx1) = [];

mov = data1.U * data1.blueV;

data2.blueV = reshape(data2.blueV, size(data2.blueV,1), []);
data2.U = arrayShrink(data2.U, mask);
trialIdx2 = isnan(data2.blueV(1, :));
data2.blueV(:, trialIdx2) = [];

mov = cat(2, mov, data2.U * data2.blueV);

disp('Computing SVD');
COV       = mov' * mov/size(mov,1);
totalVar  = sum(diag(COV)); % total variance of data.
opts.nSVD = min(size(COV,1)-2, opts.nSVD);
[V, Sv]  = rsvd(double(COV), opts.nSVD);

Sv = diag(Sv);
U = single(normc(mov * V));

blueV = NaN(opts.nSVD, length(trialIdx1) + length(trialIdx2));
blueV(:,[~trialIdx1 ~trialIdx2]) =  U' * mov;

blueV = reshape(blueV, dSize(1), dSize(2), []);
U = arrayShrink(U,mask,'split'); %recreate spatial components
save([newPath filesep 'blueV.mat'],'blueV','U','blueFrametimes','Sv','totalVar','trials','bTrials','-v7.3');

%%
opts.nSVD = 500;
fPath1 = [cPath rec1];
fPath2 = [cPath rec2];
newPath = [cPath rec3];

data1 = load([fPath1 filesep 'hemoV.mat']);
data2 = load([fPath2 filesep 'hemoV.mat']);

hemoFrameTimes = cat(2,data1.hemoFrameTimes, data2.hemoFrameTimes);
bTrials = [data1.bTrials data2.bTrials];
trials = [data1.trials data2.trials+data1.trials(end)];
mask1 = isnan(data1.U(:,:,1));
mask2 = isnan(data2.U(:,:,1));
mask = mask1 | mask2;


% kill later
data1.hemoV = data1.hemoV(:, :, 1:2);
data2.hemoV = data2.hemoV(:, :, 1:2);

dSize = size(data1.hemoV);
data1.hemoV = reshape(data1.hemoV, size(data1.hemoV,1), []);
data1.U = arrayShrink(data1.U, mask);
trialIdx1 = isnan(data1.hemoV(1, :));
data1.hemoV(:, trialIdx1) = [];

mov = data1.U * data1.hemoV;

data2.hemoV = reshape(data2.hemoV, size(data2.hemoV,1), []);
data2.U = arrayShrink(data2.U, mask);
trialIdx2 = isnan(data2.hemoV(1, :));
data2.hemoV(:, trialIdx2) = [];

mov = cat(2, mov, data2.U * data2.hemoV);

disp('Computing SVD');
COV       = mov' * mov/size(mov,1);
totalVar  = sum(diag(COV)); % total variance of data.
opts.nSVD = min(size(COV,1)-2, opts.nSVD);
[V, Sv]  = rsvd(double(COV), opts.nSVD);

Sv = diag(Sv);
U = single(normc(mov * V));

hemoV = NaN(opts.nSVD, length(trialIdx1) + length(trialIdx2));
hemoV(:,[~trialIdx1 ~trialIdx2]) =  U' * mov;

hemoV = reshape(hemoV, dSize(1), dSize(2), []);
U = arrayShrink(U,mask,'split'); %recreate spatial components
save([newPath filesep 'hemoV.mat'],'hemoV','U','hemoFrameTimes','Sv','totalVar','trials','bTrials','-v7.3');

