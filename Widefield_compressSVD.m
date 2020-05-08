function [Ur, Vr, Ulow, Vlow] = Widefield_compressSVD(data, opts, noiseTest)
% function [Ured,Vtred,meanM,keeprowids,Mlowr] = Widefield_compressSVD(M,SVD_method,maxlag,confidence,mean_threshold_factor,snr_threshold)

if ~exist('noiseTest', 'var') || isempty(noiseTest)
    noiseTest = true; %by default, there should be a test for signal vs noise dimensions
end

Ur = []; Vr = []; Ulow = []; Vlow = [];
useIdx = ~isnan(data(:,1)); %check for NaNs in the 1st dims. don't use those.
if sum(useIdx) > 0
    % run SVD
    if strcmp(opts.svdMethod, 'vanilla')
        [U, s, Vr] = svd(data, 'econ');
    elseif strcmp(opts.svdMethod, 'randomized')
        [U, s, Vr] = fsvd(data, opts.blockDims, 1, 0);
    end
    Vr = s * Vr'; %multiply S into V, so only U and V from here on
    
    if noiseTest
        % Determine which components to keep using an autocorrelations test
        ctid = choose_rank(Vr, opts.maxLag, opts.autoConfidence, opts.autoThresh);
        rankIdx = find(ctid(1, :) == 1); %dimensions to keep
        
        % Also remove those components that have too low snr
        snrIdx = (std(Vr(rankIdx,:),[],2) ./ noise_level(Vr(rankIdx,:)) > opts.snrThresh);
        rankIdx(~snrIdx) = [];
    else
        rankIdx = 1 : size(Vr,1); %if test is omitted, return all dimensions as signal dimensions
    end
    
    % split signal from noise
    selIdx = false(1, size(U,2)); selIdx(rankIdx) = true;
    Ur = NaN(size(data,1), sum(selIdx), 'single'); %preallocate reduced U that has NaNs for unused pixels
    Ulow = NaN(size(data,1), sum(~selIdx), 'single'); %preallocate reduced U that has NaNs for unused pixels

    Ur(useIdx,:) = U(useIdx, selIdx); %keep selected U
    Ulow(useIdx,:) = U(useIdx, ~selIdx); %keep rejected U
    Vlow = Vr(~selIdx, :); %keep rejected V
    Vr = Vr(selIdx, :); %keep selected V
    
    if opts.verbosity %give some feedback
        fprintf('\tInitial Rank : %d\n',sum(ctid(1, :) == 1));
        fprintf('\tlow snr components: %d \n',sum(~snrIdx));
        fprintf('\tFinal Rank : %d \n',size(Vr,1));
    end
end
end