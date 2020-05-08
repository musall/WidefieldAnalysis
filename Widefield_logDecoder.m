function [cvChoice, bMaps] = Widefield_logDecoder(data, U, leftIdx, stepSize, regType, learnType, useTrials)
%% run logistic regression decoder on widefield data
rng(1) % for reproducibility
Cnt = 0;

if isempty(U)
    bMaps = NaN;
end

for iSteps = 1 : stepSize : size(data,2)
    Cnt = Cnt + 1;
    if iSteps + stepSize - 1 < size(data,2)
        X = squeeze(nanmean(data(:, iSteps:iSteps+stepSize-1, :),2));
    else
        X = squeeze(nanmean(data(:, iSteps:end, :),2));
    end
    useIdx = ~isnan(X(1,:));  %useable trials at this time
    
    %make sure L/R responses are equally distributed
    respDiff = sum(leftIdx(useIdx)) - sum(~leftIdx(useIdx));
    if respDiff < 0 %too many rightward responses
        cIdx = find(~leftIdx & useIdx);
        cIdx = cIdx(randperm(length(cIdx), abs(respDiff))); %don't use those trials
        useIdx(cIdx) = false;
    elseif respDiff > 0  %too many leftward responses
        cIdx = find(leftIdx & useIdx);
        cIdx = cIdx(randperm(length(cIdx), abs(respDiff))); %don't use those trials
        useIdx(cIdx) = false;
    end
    
    if sum(useIdx) > useTrials
        % normalize
        X = X(:, useIdx);
        X = bsxfun(@minus,X,mean(X,2));
        Xstd = std(X,0,2);
        X = bsxfun(@rdivide,X,Xstd);
        Y = double(leftIdx(useIdx)');
        
        % get cross-validated model to compute predictive power (all trials)
        Mdl = fitrlinear(X,Y, 'ObservationsIn','columns', 'Crossval', 'on', 'Regularization',regType, 'Learner', learnType);
        cvChoice(Cnt, 1) = sum(Y == (kfoldPredict(Mdl) > 0.5)) / length(Y);
        if ~isempty(U)
            bMaps(:, Cnt, 1) = U * bsxfun(@rdivide,Mdl.Trained{1}.Beta, Xstd);
        end
    else
        cvChoice(Cnt, 1) = NaN;
        if ~isempty(U)
            bMaps(:, Cnt, 1) = NaN(size(U,1), 1);
        end
    end
end
end