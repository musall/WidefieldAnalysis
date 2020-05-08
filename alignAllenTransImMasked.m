function [im, mask] = alignAllenTransImMasked(im, transParams, fillVal)
% [im, mask] = alignAllenTransImMasked(im, transParams, fillVal)
% 
% Like alignAllenTransIm, this function rotates the image or image stack
% "im" according to the parameters in transParams (from opts2.mat). Also
% handles cleaning up the resulting image mask, using the NaNs from the
% first image in im. Masked-off values will be replaced by fillVal.


%% Prepare mask

mask = alignAllenTransIm(single(isnan(im(:, :, 1))), transParams);
mask(isnan(mask)) = 1;
mask = (mask > 0);


%% Transform images

im = alignAllenTransIm(im, transParams);

% Include out-of-rotation pixels in mask
mask = (mask | isnan(im));

% Fill-in missing values
im(mask) = fillVal;
