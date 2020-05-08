function [nV, nU, rankIdx] = Widefield_denoiseBlocks(bV, opts)

bV = cat(1,bV{:});
%merge dimensions if bV is in dims x trials x frames format
if length(size(bV)) == 3
    bV = reshape(bV,size(bV,1), []);
end

%combine all blue blocks and run a second SVD
[nU, s, nV] = fsvd(bV, min(size(bV)), 1, 0); %combine all blocks in a second SVD
nV = s * nV'; %multiply S into V, so only U and V from here on

% identify and remove noisy dimension from blue channel
temp = choose_rank(nV, opts.maxLag, opts.autoConfidence, opts.autoThresh);
rankIdx = find(temp(1, :) == 1); %dimensions to keep
bSnrIdx = (std(nV(rankIdx,:),[],2) ./ noise_level(nV(rankIdx,:))) > opts.snrThresh;
rankIdx(~bSnrIdx) = [];
nU = nU(:, rankIdx);
nV = nV(rankIdx, :);
end

function vtid=choose_rank(Vt,maxlag,confidence,mean_threshold_factor)
%%%%%%%%%%%%%%%%%%
%     Select rank wrt axcov
%%%%%%%%%%%%%%%%%%
[n, L] = size(Vt);
vtid = nan*ones(2, n);
%       Reject null hypothesis comes from white noise
mean_th = covCI_wnoise(L,confidence,maxlag);
mean_th = mean_th*mean_threshold_factor;
keep = cov_one(Vt, mean_th,maxlag);
lose = setdiff(1:n, keep);
vtid(1, keep) = 1;  % due to cov
vtid(2, lose) = 1;  % extra
end

function mean_th=covCI_wnoise(L, confidence, maxlag, n)
%%%%%%%%%%%%%%%%%%
%     Compute mean_th for auto covariance of white noise
%%%%%%%%%%%%%%%%%%
if nargin<4, n=2000; end
covs_ht = zeros(n,1);
for sample=1:n
    ht_data = randn(L,1);
    covdata = autocorr(ht_data, maxlag);
    covs_ht(sample) = mean(covdata);
end
mean_th = mean_confidence_interval(covs_ht, confidence);
end

function corr_x=autocorr(x, maxlag)
%%%%%%%%%%%%%%%%%%
%     Compute autocorrelation of vector x from autocovariance
%%%%%%%%%%%%%%%%%%
if maxlag > length(x), maxlag = length(x); end

cov_x = axcov(x, maxlag);
corr_x = cov_x/var(x);
corr_x= corr_x(maxlag+1:end);
end

function x=mean_confidence_interval(data, confidence)
%%%%%%%%%%%%%%%%%%
%     Compute mean confidence interval (CI)
%     for a normally distributed population
%     Input:
%         data: input
%         confidence: confidence level
%     Output:
%         m: mean
%         m - h: lower CI
%         m - l: upper CI
%%%%%%%%%%%%%%%%%%
pd = makedist('Normal','mu',mean(data),'sigma',std(data));
x = icdf(pd,confidence);
end

function xcov=axcov(data, maxlag)
%%%%%%%%%%%%%%%%%%
%     Compute the autocovariance of data at lag = -maxlag:0:maxlag
%
%     Parameters:
%     ----------
%     data : array
%         Array containing fluorescence data
%
%     maxlag : int
%         Number of lags to use in autocovariance calculation
%
%     Returns:
%     -------
%     axcov : array
%         Autocovariances computed from -maxlag:0:maxlag
%%%%%%%%%%%%%%%%%%

data = data - mean(data);
T = length(data);
xcov = xcorr(data,maxlag);
xcov=real(xcov)/T;
end

function keep=cov_one(Vt, mean_th, maxlag)
%%%%%%%%%%%%%%%%%%
%     Compute auto covariance, keep if mean reaches mean_th
%%%%%%%%%%%%%%%%%%
keep = []; vi_cov=zeros(size(Vt,1),maxlag+1);
for vector = 1:size(Vt,1)
    %           standarize and normalize
    vi = Vt(vector, :);
    vi = (vi - mean(vi))/std(vi);
    vi_cov(vector,:) = autocorr(vi, maxlag);
    if mean(vi_cov(vector,:)) >= mean_th
        keep=[keep vector]; %#ok<AGROW>
    end
end
end

function sns=noise_level(Vt,method,range_ff)
if nargin<2, method='logmexp'; end
if nargin<3, range_ff=[0.25,0.5]; end

sns=zeros(size(Vt,1),1);
for i = 1:size(Vt,1)
    [Pxx,ff] = pwelch(Vt(i,:),hanning(256),256/2,256,1,'onesided');
    Pxx_ind = Pxx(ff > range_ff(1) & ff < range_ff(2));
    switch method
        case 'mean'
            sns(i)=sqrt(mean(Pxx_ind/2));
        case 'median'
            sns(i)=sqrt(median(Pxx_ind/2));
        case 'logmexp'
            sns(i)=sqrt(exp(mean(log(Pxx_ind/2))));
    end
end
end