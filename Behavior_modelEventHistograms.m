function [bhvEvents, histEdges] = Behavior_modelEventHistograms(cPath,reqLabels,binSize,motorIdx,frames)
% This code is to collect behavioral events that were used as regressors in
% the linear model and reconstruct PSTHs for single trials. This is to
% assess how different kind of mouse behaviors are aligned with trial time.
% Usage: [bhvEvents, histEdges] = Behavior_modelEventHistograms(cPath,reqLabels)
%
% Inputs: cPath = path to current data folder. must contain regdata file
%                 from linear model code.
%         reqLabels = labels of requested regressors. This is a cell array
%                     of strings. Single entries can contain another cell
%                     of strings if multiple regressors should be combined
%                     into one (like lLick/rLick).
%
% Output: bhvEvents = Histograms of requested behavioral events in the
%                     order of reqlabels. 
%                     Dimensions are frames x regs x modality (vis and aud)
%         histEdges = Edges to plot histograms later.

if ~strcmpi(cPath(end),filesep)
    cPath = [cPath filesep];
end

if ~exist('binSize','var') || isempty(binSize)
    binSize = 5; % number of frames per bin
end

if ~exist('frames','var') || isempty(frames)
    frames = 205; % number of frames per bin
end

if ~exist('motorIdx','var') || isempty(motorIdx)
    motorIdx = 16; %index of zero-lag motor regressors
end

histEdges = 0:binSize:ceil(frames/binSize)*binSize; %edges for histograms

%% load data
load([cPath 'regData.mat'],'fullR', 'recIdx', 'recLabels', 'idx', 'trialIdx'); %load orthogonalized design matrix
load([cPath 'dimBeta.mat'],'dimBeta'); %load model betas for orthogonalized design matrix
load([cPath 'Vc.mat'],'bTrials'); %get index for trials that were used in the model

%% load behavior and get modality indices
bhvFile = strsplit(cPath,filesep);
bhvFile = dir([cPath bhvFile{end-3} '*' bhvFile{end-2} '*.mat']);
load([cPath bhvFile.name]);

modIdx(1,:) = reshape(repmat(SessionData.StimType(bTrials) == 1,frames,1),[],1); % visual trials
modIdx(2,:) = reshape(repmat(SessionData.StimType(bTrials) == 2,frames,1),[],1); % audio trials trials
modIdx = modIdx(:,~trialIdx); %remove non-used trials

%% get behavioral events
bhvEvents = NaN(length(histEdges)-1,length(reqLabels),2);
for x = 1 : length(reqLabels)
    vecOut = [];
    cLabel = cellstr(reqLabels{x}); %make sure input is a cell array with strings
    
    for iRegs = 1 : length(cLabel)
        cInd = find(recIdx == find(ismember(recLabels, cLabel{iRegs}))); %current regressor set
        cInd = cInd(motorIdx); %index in complete design matrix
        if idx(cInd)  %regressor got rejected in the model
            vecOut(:,iRegs) = false(size(fullR,1),1);
        else
            cInd = cInd - sum(idx(1:cInd(1))); %index in reduced design matrix
            temp = fullR(:, cInd); %current zero-mean regressor
            vecOut(:,iRegs) = logical(temp - min(temp)); %convert to logical
        end
    end
    vecOut = sum(vecOut,2); %combine if requesting a set of regressors
    
    for iMods = 1:2
        cOut = rem(find(vecOut(modIdx(iMods,:)) > 0)-1,frames); %convert to single trial 'timestamps'
        bhvEvents(:,x,iMods) = histcounts(cOut,histEdges) ./ (sum(modIdx(iMods,:)) / frames); %get histogram
    end
end
end

