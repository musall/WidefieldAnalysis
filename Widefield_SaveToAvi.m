function Widefield_SaveToAvi(DataIn, fName,fRate,cMap,Range,path,smth,sigma,imgAngle)
% Function to convert a 3D matrix into a short movie sequence. Transforms
% data matrix into an indexed image sequence that is written to .avi file.
% Usage: Widefield_SaveToAvi(DataIn, fName,fRate,cMap,Range,path,smth,sigma,imgAngle)
% Inputs: DataIn = 3D matrix of size (x,y,Frames).
%         fName = File name of .avi file (without .avi ending)
%         fRate = framerate of movie. E.g. at 10Hz, a stack of size
%                 (x,y,10) will result in a 1s long movie sequence.
%         cMap = colormap that is used. e.g. jet, gray, ect.
%         Range = Range at which colors should be used. Values below or
%                 above threshold will be equal to min(Range) or max(Range)
%                 respectively.
%       optional: path = path of .avi file. If not assigned file is written
%                        to current folder.
%                 smth = spatial smoothing. Size of a gaussian filter that is used for spatial smoothing.
%                 sigma = sigma for gaussian filter.
%                 imgAngle =  Angles in deg to rotate video if required. Rotates counter-clockwise.

if exist('path','var')==0 || isempty(path)
    path = [pwd '\'];
else
    if ~strcmp(path(end),'\')
        path = [path '\'];
    end
end
if exist('smth','var')==0 || isempty(smth)
    smth = 0;
end
if exist('sigma','var')==0 || isempty(smth)
    sigma = 0;
end
if exist('imgAngle','var') && ~isempty(imgAngle)
    DataIn = imrotate(DataIn,imgAngle);
    DataIn(~imrotate(true([size(DataIn,1) size(DataIn,2)]),imgAngle)) = NaN;
end

[x,y]=meshgrid(-smth:smth,-smth:smth); % create gaussian filter based on flength and sigma
Kernel= exp(-(x.^2+y.^2)/(2*sigma*sigma))/(2*pi*sigma*sigma);

cMap = str2func(cMap);
C = cMap(256);
C(1,:) = [0 0 0]; %make lowest value black and reserve for Nan/Inf values

v = VideoWriter([path fName]);
v.FrameRate = fRate;
% v.Quality = 100; %this can make the files a lot larger as the 75 default
open(v);
set(0,'DefaultFigureWindowStyle','docked')
h = figure;

for iFrames = 1:size(DataIn,3)
    pic = (DataIn(:,:,iFrames));
    pic = double(pic);
    if ~isnan(Kernel)
        pic = conv2(pic, Kernel,'valid');
    end
    blackInd = isnan(pic) | isinf(pic); %index for NaN/Inf values
    
    pic = mat2gray(pic,Range);
    pic(blackInd) = 0;
    pic = gray2ind(pic,256);
    pic(pic == 0 & ~blackInd) = 1;
    pic = ind2rgb(pic, C);

    writeVideo(v,pic)
end
set(0,'DefaultFigureWindowStyle','normal')
close(h);
close(v);    
    