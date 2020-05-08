function mask = Widefield_getMask(trace,iThresh)
%quick code to get a mask. Threshold based on percentile.

set(0,'DefaultFigureWindowStyle','docked')
Check = false;
h = figure('name','Check Threshold');
while ~Check
    imagesc(trace);axis square; colormap gray; hold on; drawnow
    contour(imfill(trace > prctile(trace(:),iThresh),'holes'),'w'); axis square; title(['Dark image index - Threshold: ' int2str(iThresh)])
    iThresh = prctile(trace(:),iThresh); %threshold to detect darker part of the image. Pixels below threshold are excluded from further analysis.
    Wait = input('Happy with threshold? Enter "Y" or new threshold (0-100) to proceed \n','S');
    if strcmpi(Wait,'y')
        Check = true;
    elseif  ~isempty(str2num(Wait))
        Wait = str2num(Wait);
        disp(['Changed threshold to '  num2str(Wait)]);
        iThresh = Wait;
    end
end
close(h);
mask = ~imfill(trace > iThresh, 'holes'); %mask to remove pixels
set(0,'DefaultFigureWindowStyle','normal')

end