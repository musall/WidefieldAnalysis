function mask = Widefield_CircularMask(mData)
% uses imellipse to create drawable ellipse that can be used as a mask.
h = figure('renderer','painter');
imshow(mat2gray(mData)); %create figure to draw in
h1 = imellipse;
position = wait(h1);
mask = ~poly2mask(position(:,1),position(:,2),size(mData,1),size(mData,2));
close(h);