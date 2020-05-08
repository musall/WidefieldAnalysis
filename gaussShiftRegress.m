function RMSE = gaussShiftRegress(Vc, fullR, gaussWin, frames)
% perform regression using different lengths of gaussian window (convolved 
% with design matrix) and return RMSE as indication for goodness of fit.
% Vc is data, fullR design matrix, gaussWin the tested FWHM and frames
% the number of frames in a single trial.
    
[a,b] = size(fullR);

% find non-continous regressors (contain values different from -1, 0 or 1)
temp = false(size(fullR));
temp(fullR(:) ~= 0 & fullR(:) ~= 1 & fullR(:) ~= -1 & ~isnan(fullR(:))) = true;
idx = nanmean(temp) == 0; %index for non-continous regressors

% do gaussian convolution. perform trialwise to avoid overlap across trials.
trialCnt = a/frames;
fullR = reshape(fullR,frames,trialCnt,b);
for iTrials = 1:trialCnt
    fullR(:,iTrials,idx) = smoothCol(squeeze(fullR(:,iTrials,idx)),gaussWin,'gauss');
end
fullR = reshape(fullR,a,b);

if cond(fullR) < 1000
    beta = (fullR' * fullR) \ fullR' * Vc;
    RMSE = (Vc - fullR * beta) .^ 2;
    RMSE = mean(RMSE(:));
else
    RMSE = std(Vc(:));
    warning('Design matrix is badly conditioned. Returing std(Y) for rmse.')
end
disp(RMSE);


function dataOut = smoothCol(dataIn, fWidth,fType,fLength)
% function to smooth collumns of a 2d matrix 'dataIn'.
% fType defines the filter type which is either a box or gaussian filter.
% When using a box filter 'fWidth' is the size of the moving average, when
% using a gaussian filter 'fWidth' defines its full width half maximum.
% Default filter is a 5pt box filter when only 'dataIn' is provided.
% ----------------------------------------
% To make the code run faster when using the gaussian filter, one can also
% use a shorter filter length. Default is 100*sigma.

if ~exist('fType','var')
    fType = 'box'; %default is box filter
end

if ~exist('fWidth','var')
    fWidth = 5; %default filter length is 5 dpoints
end

if ~ismember(fType,{'box' 'gauss'})
    error('Unknown filter type. Use "box" or "gauss" for fType')
end

if sum(abs(dataIn(:))) > 0 %check if there is any data available
    
    if ~strcmpi(class(dataIn),'double') %make sure data is double
        dataIn = double(dataIn);
    end
    
    if strcmpi(fType,'box') %box filter
        n = size(dataIn,1);
        cbegin = cumsum(dataIn(1:fWidth-2,:),1);
        cbegin = bsxfun(@rdivide, cbegin(1:2:end,:), (1:2:(fWidth-2))');
        cend = cumsum(dataIn(n:-1:n-fWidth+3,:),1);
        cend = bsxfun(@rdivide, cend(end:-2:1,:), (fWidth-2:-2:1)');
        dataOut = conv2(dataIn,ones(fWidth,1)/fWidth,'full'); %smooth trace with moving average
        dataOut = [cbegin;dataOut(fWidth:end-fWidth+1,:);cend];
        
    elseif strcmpi(fType,'gauss') %gaussian filter
        fSig = fWidth./(2*sqrt(2*log(2))); %in case of gaussian smooth, convert fwhm to sigma.
        if ~exist('fLength','var') || isempty(fLength)
            fLength = round(fSig * 100); %length of total gaussian filter
        end
        fLength = fLength-1+mod(fLength,2); % ensure kernel length is odd
        kernel = exp(-linspace(-fLength / 2, fLength / 2, fLength) .^ 2 / (2 * fSig ^ 2));
        dataOut = conv2(dataIn,kernel','same'); %smooth trace with gaussian
    end
    
else
    dataOut = dataIn;
end