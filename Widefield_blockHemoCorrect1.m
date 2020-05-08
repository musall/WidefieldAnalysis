function [sNewU, sNewV, avgNewBlock] = Widefield_blockHemoCorrect1(sBlueV, sBlueU, sHemoV, sHemoU, blueFrameTimes, opts)
% Hemodynamic correction after running blockwise SVD and denoising. This
% code assumes that blue and hemo channel have different V/Us and computes
% a new, corrected V/U by reconstructing each block and running another
% denoiser step.

smth = ceil(0.15 / (1/opts.sRate)); %use 175ms window for moving average when smoothing hemoV

for iBlocks = 1 : size(sBlueV,1)
    
    
    for iTrials = 1 : length(blueFrameTimes)
        temp = hemoBlock(:, Cnt + (1 : length(blueFrameTimes{iTrials})));
        
        cbegin = cumsum(temp(:,1:smth-2),2);
        cbegin = bsxfun(@rdivide, cbegin(:,1:2:end), 1:2:(smth-2));
        
        cend = cumsum(temp(:,end:-1:end-smth+3),2);
        cend = bsxfun(@rdivide, cend(:,end:-2:1), (smth-2:-2:1));
        
        temp = conv2(temp,ones(1,smth)/smth,'full'); %smooth trace with moving average of 'smth' points
        hemoBlock(:, Cnt + (1 : length(blueFrameTimes{iTrials}))) = [cbegin temp(:,smth:end-smth+1) cend];
        Cnt = Cnt + length(blueFrameTimes{iTrials});
    end
    

    
    
    
    
    
    % reconstruct current block and high-pass filter above 0.1Hz
    blueBlock = sBlueU{iBlocks} * sBlueV{iBlocks};
    hemoBlock = sHemoU{iBlocks} * sHemoV{iBlocks};
    
    [b, a] = butter(2,0.1/(opts.sRate), 'high');
    blueBlock = single(filtfilt(b,a,double(blueBlock')))';
    hemoBlock = single(filtfilt(b,a,double(hemoBlock')))';
   
    % smooth hemo V
    Cnt = 0;
    smth = smth - 1 + mod(smth,2); % ensure kernel length is odd
    for iTrials = 1 : length(blueFrameTimes)
        temp = hemoBlock(:, Cnt + (1 : length(blueFrameTimes{iTrials})));
        
        cbegin = cumsum(temp(:,1:smth-2),2);
        cbegin = bsxfun(@rdivide, cbegin(:,1:2:end), 1:2:(smth-2));
        
        cend = cumsum(temp(:,end:-1:end-smth+3),2);
        cend = bsxfun(@rdivide, cend(:,end:-2:1), (smth-2:-2:1));
        
        temp = conv2(temp,ones(1,smth)/smth,'full'); %smooth trace with moving average of 'smth' points
        hemoBlock(:, Cnt + (1 : length(blueFrameTimes{iTrials}))) = [cbegin temp(:,smth:end-smth+1) cend];
        Cnt = Cnt + length(blueFrameTimes{iTrials});
    end
    
    % subtract mean of all pixels and compute single-pixel regression coefficients.
    blueBlock = bsxfun(@minus, blueBlock, mean(blueBlock,2));
    hemoBlock = bsxfun(@minus, hemoBlock, mean(hemoBlock,2));
    regC = nansum(blueBlock.*hemoBlock,2) ./ nansum(hemoBlock.*hemoBlock,2);

    % compute new data and compress again
    newBlock = blueBlock - (hemoBlock .* regC);
    avgNewBlock{iBlocks} = mean(newBlock, 2);    
    [sNewU{iBlocks}, sNewV{iBlocks}] = Widefield_compressSVD(newBlock, opts, false);
    
end