function im = alignAllenTransIm(im, transParams)
% im = alignAllenTransIm(im, transParams)
% 
% Take a brain image and rotate, scale, and translate it according to the
% parameters in transParams (produced by alignBrainToAllen GUI). im may be
% an image stack. Note that due to the behavior of imrotate, NaNs may not
% propagate identically if an image stack is used vs. calling this function
% on individual images.
% 
% Pixels where the value is not defined (due to rotation) are set to NaN.
% This behavior is different from imrotate.

offset = 5E1;

dSize = size(im);
[h, w, d] = size(im);

% Set pixels off from zero
[theMin, minIdx] = min(im(:));
im = im - theMin + offset;

% Rotate
im = imrotate(im, transParams.angleD, 'bilinear');

% Scale
if transParams.scaleConst ~= 1
  im = imresize(im, transParams.scaleConst);
end

% Set NaNs to 0, because imtranslate can't handle them
nans = isnan(im);
if any(nans(:))
  im(nans) = 0;
end

% Translate
im = imtranslate(im, transParams.tC');

% Detect missing pixels due to rotation, set to NaN
im(im < offset) = NaN;

% Restore offset
im = im + theMin - offset;
im(minIdx) = theMin;

% Trim result (rotate expands it)
newSize = size(im);
trimH = floor((newSize(1) - dSize(1)) / 2);
trimW = floor((newSize(2) - dSize(2)) / 2);
im = im(trimH + (1:dSize(1)), trimW + (1:dSize(2)), :);
im = reshape(im, dSize);
