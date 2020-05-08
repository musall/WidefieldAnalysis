function [cV, regC, T, hemoExp] = Widefield_blockHemoCorrect(bV, bU, blueFrameTimes, opts)
% Hemodynamic correction after running blockwise SVD and denoising. This
% code assumes that blue and hemo channel have different V/Us and computes
% a new, corrected V/U by reconstructing each block and running another
% denoiser step.

smth = ceil(0.15 / (1/opts.sRate)); %use 175ms window for moving average when smoothing hemoV
[b, a] = butter(2,1/2, 'low');

for iBlocks = 1 : size(bV,1)
    
    % split into two channels
    cBlock = reshape(bV{iBlocks}, size(bV{iBlocks},1), [], 2);
    blueV = cBlock(:, :, 1);
    hemoV = cBlock(:, :, 2);
    clear cBlock

    % smooth hemo V
    Cnt = 0;
    smth = smth - 1 + mod(smth,2); % ensure kernel length is odd
    for iTrials = 1 : length(blueFrameTimes)
        temp = blueV(:, Cnt + (1 : length(blueFrameTimes{iTrials})))';
        temp = [repmat(temp(1,:),10,1); temp; repmat(temp(end,:),10,1)];
        temp = single(filtfilt(b,a,double(temp)))';
        blueV(:, Cnt + (1 : length(blueFrameTimes{iTrials}))) = temp(:, 11:end-10);
        
        temp = hemoV(:, Cnt + (1 : length(blueFrameTimes{iTrials})))';
        temp = [repmat(temp(1,:),10,1); temp; repmat(temp(end,:),10,1)];
        temp = single(filtfilt(b,a,double(temp)))';
        hemoV(:, Cnt + (1 : length(blueFrameTimes{iTrials}))) = temp(:, 11:end-10);
        Cnt = Cnt + length(blueFrameTimes{iTrials});
    end
    
    % reconstruct current blocks
    blueBlock = bU{iBlocks} * blueV;
    hemoBlock = bU{iBlocks} * hemoV;
    
    % subtract mean of all pixels and compute single-pixel regression coefficients.
    blueBlock = bsxfun(@minus, blueBlock, mean(blueBlock,2));
    hemoBlock = bsxfun(@minus, hemoBlock, mean(hemoBlock,2));
    regC = nansum(blueBlock.*hemoBlock,2) ./ nansum(hemoBlock.*hemoBlock,2);
    
    % make the prediction
    T = pinv(bU{iBlocks}) * bsxfun(@times, regC(:), bU{iBlocks});
    cV{iBlocks} = (blueV' - hemoV'*T')';
    
    % compute variance explained
    bVar(iBlocks) = nansum(blueBlock(:).^2);
    cVar(iBlocks) = nansum(cV{iBlocks}(:).^2);

end
hemoExp = 100*(sum(bVar)-sum(cVar))/sum(bVar);
fprintf('%f percent variance explained by hemodynamic correction\n', hemoExp);
    
end