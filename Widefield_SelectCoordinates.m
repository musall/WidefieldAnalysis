function [x, y] = Widefield_SelectCoordinates(dataIn,pxPerMM,target)
        
grid(gca,'on');grid minor;set(gca,'GridColor','w');
Check = false;
while ~Check
    title(['Select bregma ' target ' to align grid']);
    [x,y] = ginput(1);
    xVec = x-floor(x/pxPerMM)*pxPerMM:pxPerMM:size(dataIn,1);
    xLabel = num2str(abs((1:length(xVec))-ceil(x/pxPerMM))');
    set(gca,'xTick',xVec); set(gca,'xTickLabel',xLabel)
    
    yVec = y-floor(y/pxPerMM)*pxPerMM:pxPerMM:size(dataIn,2);
    yLabel = num2str(abs((1:length(yVec))-ceil(y/pxPerMM))');
    set(gca,'yTick',yVec); set(gca,'yTickLabel',yLabel)
    
    Wait = input(['Type "y" to continue or any other key to set ' target ' again. \n '],'S');
    if strcmpi(Wait,'y')
        Check = true;
    end
end