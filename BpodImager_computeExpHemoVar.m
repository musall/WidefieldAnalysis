function BpodImager_computeExpHemoVar
sRate = 30;
figure
%% GFP data
[dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
controlType = 'GFP';
dataOverview = dataOverview(ismember(dataOverview(:,4), controlType), :);

for iAnimal = 1 : size(dataOverview,1)
    
opts.fPath = [cPath dataOverview{iAnimal,1} filesep 'SpatialDisc' filesep dataOverview{iAnimal,3} filesep];
load([opts.fPath 'blueV'], 'blueV');
load([opts.fPath 'Vc'], 'Vc');
blueV = blueV(1:200,:,:);
Vc = Vc(1:200,:,:);

[A,B,C] = size(blueV);
blueV = reshape(blueV,A,[])';

% subtract means
blueV = bsxfun(@minus, blueV, nanmean(blueV));

% high-pass blueV and hemoV above 0.1Hz
[b, a] = butter(2,0.1/(sRate), 'high');
blueV(~isnan(blueV(:,1)),:) = single(filtfilt(b,a,double(blueV(~isnan(blueV(:,1)),:))));

% smooth hemo V
[b, a] = butter(2,15/(sRate), 'low');
blueV = reshape(blueV',A,B,C);

for iTrials = 1:C
    cIdx = ~isnan(blueV(1,:,iTrials)); %make sure to only use non-NaN frames
    blueV(:,cIdx,iTrials) = single(filtfilt(b,a,double(blueV(:,cIdx,iTrials)')))';
end

f1Pow = nansum(blueV(:).^2);
f1Powcor = nansum(Vc(:).^2);
GFP_HemoVar(iAnimal,1) = f1Pow;
GFP_HemoVar(iAnimal,2) = f1Powcor;
GFP_HemoVar(iAnimal,3) = 100*(f1Pow-f1Powcor)/f1Pow;
disp(opts.fPath);
fprintf('%f percent variance explained by hemo signal\n', 100*(f1Pow-f1Powcor)/f1Pow);
end

plot(GFP_HemoVar(:,3), 'bo-', 'linewidth', 2); hold on;

%% camK data
[dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings_GFP;
controlType = 'gcamp';
dataOverview = dataOverview(ismember(dataOverview(:,4), controlType), :);

for iAnimal = 1 : size(dataOverview,1)
    
opts.fPath = [cPath dataOverview{iAnimal,1} filesep 'SpatialDisc' filesep dataOverview{iAnimal,3} filesep];
load([opts.fPath 'blueV'], 'blueV');
load([opts.fPath 'Vc'], 'Vc');
blueV = blueV(1:200,:,:);
Vc = Vc(1:200,:,:);

[A,B,C] = size(blueV);
blueV = reshape(blueV,A,[])';

% subtract means
blueV = bsxfun(@minus, blueV, nanmean(blueV));

% high-pass blueV and hemoV above 0.1Hz
[b, a] = butter(2,0.1/(sRate), 'high');
blueV(~isnan(blueV(:,1)),:) = single(filtfilt(b,a,double(blueV(~isnan(blueV(:,1)),:))));

% smooth hemo V
[b, a] = butter(2,15/(sRate), 'low');
blueV = reshape(blueV',A,B,C);

for iTrials = 1:C
    cIdx = ~isnan(blueV(1,:,iTrials)); %make sure to only use non-NaN frames
    blueV(:,cIdx,iTrials) = single(filtfilt(b,a,double(blueV(:,cIdx,iTrials)')))';
end

f1Pow = nansum(blueV(:).^2);
f1Powcor = nansum(Vc(:).^2);
camK_HemoVar(iAnimal-5,1) = f1Pow;
camK_HemoVar(iAnimal-5,2) = f1Powcor;
camK_HemoVar(iAnimal-5,3) = 100*(f1Pow-f1Powcor)/f1Pow;
disp(opts.fPath);
fprintf('%f percent variance explained by hemo signal\n', 100*(f1Pow-f1Powcor)/f1Pow);
end

plot(camK_HemoVar(:,3), 'ro-', 'linewidth', 2); hold on;

%% WF data
[dataOverview, ~, ~, ~, ~, ~, ~, cPath] = delayDecRecordings;
% cPath = 'X:\smusall\BpodImager\Animals\';

for iAnimal = 1 : size(dataOverview,1)
    
opts.fPath = [cPath dataOverview{iAnimal,1} filesep 'SpatialDisc' filesep dataOverview{iAnimal,3} filesep];
cd(opts.fPath);
load('blueV.mat', 'blueV');
load('Vc.mat', 'Vc');
blueV = blueV(1:200,:,:);
Vc = Vc(1:200,:,:);

[A,B,C] = size(blueV);
blueV = reshape(blueV,A,[])';

% subtract means
blueV = bsxfun(@minus, blueV, nanmean(blueV));

% high-pass blueV and hemoV above 0.1Hz
[b, a] = butter(2,0.1/(sRate), 'high');
blueV(~isnan(blueV(:,1)),:) = single(filtfilt(b,a,double(blueV(~isnan(blueV(:,1)),:))));

% smooth hemo V
[b, a] = butter(2,15/(sRate), 'low');
blueV = reshape(blueV',A,B,C);

for iTrials = 1:C
    cIdx = ~isnan(blueV(1,:,iTrials)); %make sure to only use non-NaN frames
    blueV(:,cIdx,iTrials) = single(filtfilt(b,a,double(blueV(:,cIdx,iTrials)')))';
end

f1Pow = nansum(blueV(:).^2);
f1Powcor = nansum(Vc(:).^2);
WF_HemoVar(iAnimal,1) = f1Pow;
WF_HemoVar(iAnimal,2) = f1Powcor;
WF_HemoVar(iAnimal,3) = 100*(f1Pow-f1Powcor)/f1Pow;
disp(opts.fPath);
fprintf('%f percent variance explained by hemo signal\n', 100*(f1Pow-f1Powcor)/f1Pow);
end

plot(WF_HemoVar(:,3), 'ko-', 'linewidth', 2); hold on;
