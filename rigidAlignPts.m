function [R, t] = rigidAlignPts(P, Q, allowScaling, W)
% [R, t] = rigidAlignPts(P, Q [, allowScaling] [, W])
% 
% Compute the optimal rotation and translation to align the points in P to
% the points in Q. Points are columns. P and Q must have the same number of
% points, so they should be the same size. "Optimal" is in the
% least-squares sense.
% 
% Optionally, rescaling the data may be allowed by setting allowScaling to
% 1. In this case, R includes the scaling. Default 0.
% 
% Optionally, points may be weighted by supplying a vector W of length
% size(P, 2). These values should be real, non-NaN, and >=0. If any NaNs
% are included in W, W will be ignored.
% 
% The result should have Q == R * P + t
% 
% This is an implementation of Umeyama's generalization of the Kabsch
% algorithm.


%% Error checking

if size(P) ~= size(Q)
  error('rigidAlignPts:sizeMismatch', 'P and Q must be the same size');
end


%% Optional arguments

if ~exist('allowScaling', 'var')
  allowScaling = 0;
end

if exist('W', 'var') && ~any(isnan(W))
  % More error checking
  if ~isvector(W) || length(W) ~= size(P, 2)
    error('rigidAlignPts:badWLength', ...
      'Weights vector must contain as many weights as points in P and Q');
  end
  
  if any(W < 0 | ~isreal(W))
    error('rigidAlignPts:badWSign', ...
      'Weights vector must not contain real values > 0');
  end
  
  % Ensure W is a row vector
  W = W(:)';
  
  weighted = 1;
else
  weighted = 0;
end


%% Center the data

if weighted
  meanP = bsxfun(@times, P, W) / sum(W);
  meanQ = bsxfun(@times, Q, W) / sum(W);
else
  meanP = mean(P, 2);
  meanQ = mean(Q, 2);
end

X = bsxfun(@minus, P, meanP);
Y = bsxfun(@minus, Q, meanQ);


%% Compute the rotation matrix

% Cov matrix
if weighted
  C = X * diag(W) * Y';
else
  C = X * Y';
end

[U, S, V] = svd(C);

% Because of possible mirroring in SVD, need to use determinant to orient
% correctly
orient = ones(1, size(C, 1));
orient(end) = sign(det(V * U'));
orient = diag(orient);

R = V * orient * U';


%% Compute scale

if allowScaling
  sigma2 = sum(X(:) .^ 2);
  scale = (1 ./ sigma2) * trace(S * orient);
  R = scale * R;
end


%% Compute translation vector

t = meanQ - R * meanP;
