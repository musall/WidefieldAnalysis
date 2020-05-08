function [aInd,mask] = Widefield_DrawOutline(DataIn, pThresh, mode, DoPlot, aColor, Range)

%% check variables
if ~exist('pThresh','var')
    pThresh = 5; % trehold for area detection
end
if ~exist('mode','var')
    mode = 'auto'; % mode of threshold detection
end
if ~exist('aColor','var')
    aColor = 'k'; %default color for plotline
end
if ~exist('Range','var')
    Range = []; % range for colorspace
end
if ~exist('DoPlot','var')
    DoPlot = false; %flag to plot area outline
end
tData = DataIn; %for later ploting

%% do things
if strcmpi(mode, 'manual') ||  strcmpi(mode, 'subregion') %manual selection or pre-select a subregion of the image for automatic thresholding
    
    h = figure('name','select subregion');
    DataIn(isinf(DataIn)) = 0;
    imagesc(DataIn);colormap('jet'); %create figure to draw in
    axis square; colorbar;
    if ~isempty(Range)
        caxis([Range(1) Range(2)]);
    end
    
    %% initiate drawing
    [xPos,yPos,button] = ginput(1);
    hold on
    LinePlot = plot(xPos(1),yPos(1),'b','linewidth',2);
    
    count=1;
    while button~=3 && button~=27
        count=count+1;
        [xPos(count),yPos(count),button]=ginput(1);
        if button~=3 && button~=27
            set(LinePlot,'xdata',xPos,'ydata',yPos,'LineWidth',2,'Color',aColor); %draw line to new datapoint
        else %if right mouse button or esc is pressed
            xPos(count)=[];yPos(count)=[]; %discard datapoint if right mouse or esc was pressed
        end
    end
    
    %% create index vector
    [x,y] = size(DataIn);
    xInd = ones(x,1)*(1:y);yInd = (1:x)'*ones(1,y);
    mask = inpolygon(xInd,yInd,xPos,yPos); %logical mask for ROI
    close(h);
    
    %% check if code should proceed to auto detection
    if strcmpi(mode, 'subregion')
        mode = 'auto';
    end
end

if strcmpi(mode,'auto')
    %% find areas
    C = []; Thresh = inf;
    while isempty(C)
        DataIn = tData;
        tInd = DataIn>Thresh;
        DataIn(tInd) = NaN;
        tData = DataIn;
        DataIn(~mask) = NaN; %ignore data that is outside of manually selected mask
        Thresh = prctile(reshape(DataIn,1,numel(DataIn)),pThresh); %find absolute threshold
        tInd = DataIn>Thresh; 
        C = contourc(double(tInd),1);
    end
        
    xdata = C(1,:);   %not really useful, in most cases delimters are not clear
    ydata = C(2,:);   %therefore further steps to determine the actual curves
    
    %% find curves
    n(1) = 1;         %indices where the certain curves start
    d(1) = ydata(1);  %distance to the next index
    ii = 1;
    while true
        n(ii+1) = n(ii)+d(ii)+1;    %calculate index of next startpoint
        if n(ii+1) > numel(xdata)   %breaking condition
            n(end) = [];            %delete breaking point
            break
        end
        d(ii+1) = ydata(n(ii+1));   %get next distance
        ii = ii+1;
    end
    
    [~,b]=max(d); %startpoint of largest area
    if b == length(n)
        aInd = n(b)+1:length(xdata); %index for x and y data for area if last index is used
    else
        aInd = n(b)+1:n(b+1)-1; %index for x and y data for area
    end

    %% create index vector
    [y,x,~] = size(DataIn);
    xInd = ones(y,1)*(1:x); yInd = (1:y)'*ones(1,x);
    mask = inpolygon(xInd,yInd,xdata(aInd),ydata(aInd)); %logical mask for ROI
end
aInd = contourc(double(mask),1); %index vector to draw polygon
aInd(:,aInd(1,:) < 1) = [];aInd(:,aInd(2,:) < 1) = [];
    
%% Plot area
if DoPlot
    %% select frames of interest and check for selection
    figure('name',['Selection type: ' mode]);
    imagesc(tData);hold on;
    plot(aInd(1,:),aInd(2,:),aColor,'linewidth',3);colorbar;colormap('jet');
    axis square;
    if ~isempty(Range)
        caxis([Range(1) Range(2)]);
    end
end