function [mask,outline,areaProps] = getAreas(DataIn,thresh,type,threshSign)
%function to identiy a given area that is identified by thresholding.
%DataIn is a 2D matrix, thresh is the threshold in standard deviations that is used to
%find areas and type is the property that is used to select an area. Type
%'help regionprops' to identify properties that can be used. Standard is
%'Area' which identifies the area size and will return the largest area.

if ~exist('type','var')
    type = 'Area';  % Property that is used to select area. Area is size of the area in pixels.
end

if ~exist('threshSign','var')
    threshSign = true;  % Sign of threshold. If true, the code finds an area above a given threshold.
end

if threshSign
    ind = imfill(DataIn > nanstd(DataIn(:))*thresh, 'holes'); %threshold input matrix and make sure there are no holes in separate areas
else
    ind = imfill(DataIn < nanstd(DataIn(:))*thresh, 'holes'); %threshold input matrix. if threshSign is false, search for area below the threshold.
end

foundAreas = bwlabel(ind, 8); %find independent areas
areaProps = regionprops(ind, DataIn,'all'); %get area properties

for iArea = 1:size(areaProps,1)
    sortProp(iArea) = areaProps(iArea).(type);
end

areaOutline = bwboundaries(ind);
[~,ind] = max(sortProp); %find area that maximizes the input property (if type = 'Area', this will select the largest area).
mask = foundAreas == ind;
outline = areaOutline{ind};
