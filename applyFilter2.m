function DataOut = applyFilter2(DataIn,fLength,sigma,type)
% Apply either gaussian or box filter. 
% Usage: DataOut = ApplyGauss(DataIn,fLength,sigma,type)
% Input:    DataIn = 2D matrix on which the filter should be applied
%           fLength = filter length.
%           sigma = sigma of the gaussian kernel. Standard is 1.76 for a single standard deviation.
%           type: type of filter. 1 for gaussian, 2 for box filter
% Output:   DataOut = filtered version of DataIn.

if fLength > 1
    if type == 1
        Kernel = ones(fLength,fLength); % Create box filter based on flength
        Kernel = Kernel ./ fLength^2;
    elseif type == 2
        [x,y]=meshgrid(-fLength:fLength,-fLength:fLength); % create gaussian filter based on flength and sigma
        Kernel= exp(-(x.^2+y.^2)/(2*sigma*sigma))/(2*pi*sigma*sigma);
    end
    DataOut = conv2(double(DataIn), Kernel, 'same'); %convolve original matrix with filter
else
    DataOut = DataIn; %don't apply filter if fLength is equal or lower then 1
end
% DataOut = padarray(DataOut,[1 1],NaN,'both'); %add NaNs to image to see its edge against zero-padded background when rotating.