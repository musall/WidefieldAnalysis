function mask = Widefield_ManualMask(mData)
aColor = 'w';
h = figure('name','manual selection', 'renderer', 'painter');
imshow(mData); %create figure to draw in
caxis([0 prctile(mData(:),90)*.75]);

%% initiate drawing
[xPos,yPos,button] = ginput(1);
hold on
LinePlot = plot(xPos(1),yPos(1),'b','linewidth',1);

count=1;
while button~=3 && button~=27
    count=count+1;
    [xPos(count),yPos(count),button]=ginput(1);
    if button~=3 && button~=27
        set(LinePlot,'xdata',xPos,'ydata',yPos,'LineWidth',1,'Color','w'); %draw line to new datapoint
    else %if right mouse button or esc is pressed
        xPos(count)=[];yPos(count)=[]; %discard datapoint if right mouse or esc was pressed
    end
end

%% create index vector and close figure
[x,y] = size(mData);
xInd = ones(x,1)*(1:y);yInd = (1:x)'*ones(1,y);
mask = inpolygon(xInd,yInd,xPos,yPos); %logical mask for ROI
close(h);