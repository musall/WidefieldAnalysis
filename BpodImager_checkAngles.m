
cPath = 'U:\space_managed_data\BpodImager\Animals\';
dataOverview = delayDecRecordings;
animals = dataOverview(:,1);
recs = dataOverview(:,3);
figure; hold on
for iAnimals = 1 : length(animals)
    tic
    fPath = [cPath animals{iAnimals} filesep 'SpatialDisc' filesep recs{iAnimals} filesep];
    disp(fPath); %current data path
    load([fPath 'regData.mat'], 'fullR');
    
    [~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize design matrix
    plot(abs(diag(fullQRR))); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
    if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
        error('Design matrix is rank-defficient')
    end
    
end
axis square