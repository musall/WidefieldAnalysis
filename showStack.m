function varargout = showStack(varargin)
% Function to visualize a 3D data stack. 
% Usage: showStack(matrix);
%        Input: A 3-D matrix of arrangement X,Y,Frames to scroll through.  
%        If input is a 4-D matrix,'showstack' will average over the last dimension.
%
%       GUI inputs:
%       Colorrange is given by min and max color. 
%       The slider at the figure bottom is used to scroll through frames in the stack. 
%       Filter length determines the size of a filter that can be applied
%       to the image. The filter type is either a gaussian or a box filter
%       and can be chosen with the drop-out menu at the bottom left. 
%       If chosen filter is gaussian, the sigma of the kernel is determined
%       by 'Filter sigma'.
%       Smoothed frames can be saved to base workspace using the 'Save to
%       workspace' button in the top left button. Depending on filter
%       length, the smoothed frame will be slightly smaller as the original
%       frame because the code only uses the 'valid' part of the frame
%       without zero-padded edges.
%

% If you experience unexpected behavior of the code or have suggestions for improvment, 
% leave a comment on FEX
%
% showstack was created, using GUIDE in Matlab2014b. 
% sm 5/8/2016

%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showStack_OpeningFcn, ...
                   'gui_OutputFcn',  @showStack_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before showStack is made visible.
function showStack_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to showStack (see VARARGIN)

% Choose default command line output for showStack
handles.output = hObject;

if isempty(varargin)
    error('No input data found. Data should be a 3- or 4-D frame stack')
end

handles.UserData.Close = false;
if length(size(varargin{1})) == 4
    handles.UserData.Frames = gather(squeeze(nanmean(varargin{1},4))); %4D array. Should be X,Y,Frames - Average over last dimension.
elseif length(size(varargin{1})) == 3 || length(size(varargin{1})) == 2
    handles.UserData.Frames = gather(varargin{1}); %3D / 2D array. Should be X,Y,Frames/Trials
else
    error('Wrong input. Data should be a 2, 3- or 4-D frame stack')
end

%initialize frame image
fLength = round(str2double(get(handles.FrameSmth,'string'))); %determines degree of smoothing.
data = ApplyFilter2(squeeze(handles.UserData.Frames(:,:,1)),fLength,str2double(handles.FrameSigma.String),get(handles.FilterType,'value')); %get first frame to display
rotAngle = str2double(handles.rotAngel.String); %angle of rotation in degrees. 
if ~isinf(rotAngle) % Rotate image
    data = imrotate(data,rotAngle,'crop');
    mask = ~imrotate(true(size(data)),rotAngle,'crop');
    data(mask) = NaN;
end
cImg = imshow(data,[str2double(get(handles.MinColor,'String')) str2double(get(handles.MaxColor,'String'))],'parent',handles.FrameImage); %create image object for preview
set(cImg,'AlphaData',~isnan(data)); %make NaNs transparent.
axis(handles.FrameImage,'image');
set(handles.CurrentFrame,'string',['Showing Frame 1 / ' num2str(size(handles.UserData.Frames,3))])

if handles.UseThreshold.Value
    ind = thresholdImage(data,str2double(handles.ImgThresh.String));
    set(cImg,'AlphaData',ind); %make values below transparent. Since there is no image below, pixels should appear white.
end
handles.FrameImage.Tag = 'FrameImage'; %add tag so axes can be identified later

if handles.UseThreshold.Value
    ind = thresholdImage(data,str2double(handles.ImgThresh.String));
    set(cImg,'AlphaData',ind); %make values below transparent. Since there is no image below, pixels should appear white.
end

set(handles.FrameSlider,'Min',1);
set(handles.FrameSlider,'Max',size(handles.UserData.Frames,3));
set(handles.FrameSlider,'Value',1);

if size(handles.UserData.Frames,3) == 1
    set(handles.FrameSlider,'SliderStep',[0 0])
else
    set(handles.FrameSlider,'SliderStep',[1/(size(handles.UserData.Frames,3)-1) 1/(size(handles.UserData.Frames,3)-1)]);
end

cRange = [min(data(:)) max(data(:))];
if ~any(isnan(cRange))
    caxis(handles.FrameImage,[min(data(:)) max(data(:))]);
    handles.MinColor.String = num2str(round(double(min(data(:))),4));
    handles.MaxColor.String = num2str(round(double(max(data(:))),4));
end

colormap(handles.FrameImage,strtrim(handles.pickColormap.String{handles.pickColormap.Value}));

handles.areaProp.String = fieldnames(regionprops(1,'all')); %get area properties
handles.areaProp.Value = 1; %set 1 as default (should be 'area')


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes showStack wait for user response (see UIRESUME)
uiwait(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = showStack_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
[varargout{1},varargout{2},varargout{3},varargout{4},varargout{5}] = SaveFrame_Callback([],[],handles);

if handles.UserData.Close
    delete(handles.ShowStack)
else
    handles.UserData.Close = true; %allow close request function to close figure
    guidata(hObject, handles);
end

% --- Executes on slider movement.
function FrameSlider_Callback(hObject, eventdata, handles)
% hObject    handle to FrameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(hObject,'Value',round(get(hObject,'Value')));
set(handles.CurrentFrame,'string',['Showing Frame ' num2str(get(hObject,'Value')) ' / ' num2str(size(handles.UserData.Frames,3))])
drawImage(handles,handles.FrameImage,handles.FrameSlider);


% --- Executes during object creation, after setting all properties.
function FrameSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function MaxColor_Callback(hObject, eventdata, handles)
% hObject    handle to MaxColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxColor as text
%        str2double(get(hObject,'String')) returns contents of MaxColor as a double

minColor = str2double(handles.MinColor.String);
maxColor = str2double(handles.MaxColor.String);

if isnan(maxColor) || minColor > maxColor
   handles.MaxColor.String = num2str(minColor + 1);
   maxColor = minColor + 1;
end

    
% --- Executes during object creation, after setting all properties.
function MaxColor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MinColor_Callback(hObject, eventdata, handles)
% hObject    handle to MinColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinColor as text
%        str2double(get(hObject,'String')) returns contents of MinColor as a double

minColor = str2double(handles.MinColor.String);
maxColor = str2double(handles.MaxColor.String);

if isnan(minColor) || minColor > maxColor
   handles.MinColor.String = num2str(maxColor - 1);
   minColor = maxColor -1;
end


% --- Executes during object creation, after setting all properties.
function MinColor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function FrameSmth_Callback(hObject, eventdata, handles)
% hObject    handle to FrameSmth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameSmth as text
%        str2double(get(hObject,'String')) returns contents of FrameSmth as a double

if ~isnan(str2double(get(hObject,'string')))
    FilterType_Callback([],[],handles);
else
    set(hObject,'string','0')
end

% --- Executes during object creation, after setting all properties.
function FrameSmth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameSmth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in FilterType.
function FilterType_Callback(hObject, eventdata, handles)
% hObject    handle to FilterType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FilterType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FilterType
   
drawImage(handles,handles.FrameImage,handles.FrameSlider);


% --- Executes during object creation, after setting all properties.
function FilterType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FilterType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function FrameSigma_Callback(hObject, eventdata, handles)
% hObject    handle to FrameSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FrameSigma as text
%        str2double(get(hObject,'String')) returns contents of FrameSigma as a double

if isnan(str2double(get(hObject,'string')))
    set(hObject,'string','1.76');
    disp([get(hObject,'string') ' is not a valid input for sigma. Set back to default (1.76)'])
else
    FilterType_Callback([],[],handles);
end

% --- Executes during object creation, after setting all properties.
function FrameSigma_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FrameSigma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in SaveFrame.
function [ImageOut,areaBounds,areaMask,traceData,tracePosition] = SaveFrame_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cFrame = ApplyFilter2(squeeze(handles.UserData.Frames(:,:,get(handles.FrameSlider,'Value'))),round(str2double(get(handles.FrameSmth,'string'))),str2double(handles.FrameSigma.String),get(handles.FilterType,'value')); %smooth current frame
rotAngle = str2double(handles.rotAngel.String); %angle of rotation in degrees.
if ~isinf(rotAngle) % Rotate image
    cFrame = imrotate(cFrame,rotAngle,'crop');
    mask = ~imrotate(true(size(cFrame)),rotAngle,'crop');
    cFrame(mask) = NaN;    
end
ImageOut = cFrame;
assignin('base','ImageOut',ImageOut); %send current image of both axes in cell array. Cell 1 is leftaxis , cell 2 is right axis.

lines = findall(handles.FrameImage.Children,'Type','line'); %send area outlines to base workspace
areaBounds = {}; areaMask = [];
for x = 1:length(lines)
    areaBounds{x}(1,:) = lines(x).XData;
    areaBounds{x}(2,:) = lines(x).YData;
    areaMask{x} = poly2mask(areaBounds{x}(1,:),areaBounds{x}(2,:),size(cFrame,1),size(cFrame,2));
end
if ~isempty(areaBounds)
    assignin('base','areaBounds',areaBounds);
    assignin('base','areaMask',areaMask);
end

traceOutlines = findall(handles.FrameImage.Children,'Type','rectangle'); %send area outline to base workspace
tracePosition = [];traceData = [];

if ishandle(99) %check if figure with frame traces is open.
    traces = findall(99,'Type','line');
else
    traces = [];
end

for x = 1:length(traceOutlines)
    tracePosition(x,:) = traceOutlines(x).Position; %position of rectangle that was used to create trace
    if length(traces) >= x
        traceData(x,:) = traces(x).YData; %trace data. changes in signal strenght across frames.
    end
end
if ~isempty(traceData)
    assignin('base','traceData',traceData);
end
if ~isempty(tracePosition)
    assignin('base','tracePosition',tracePosition);
end

% --- Executes on button press in SaveFrameBuffer.
function SaveFrameBuffer_Callback(hObject, eventdata, handles)
% hObject    handle to SaveFrameBuffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function pointerSize_Callback(hObject, eventdata, handles)
% hObject    handle to pointerSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pointerSize as text
%        str2double(get(hObject,'String')) returns contents of pointerSize as a double


% --- Executes during object creation, after setting all properties.
function pointerSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pointerSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in getTrace.
function getTrace_Callback(hObject, eventdata, handles)
% hObject    handle to getTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if size(handles.UserData.Frames,3) > 1
    datacursormode on
    cursorMode = datacursormode(handles.ShowStack);
    currentCursor = getCursorInfo(cursorMode);
    while isempty(currentCursor)
        drawnow;
        currentCursor = getCursorInfo(cursorMode);
        pause(0.01);
    end
    
    if strcmp(currentCursor.Target.Parent.Tag,'FrameImage')
    cPos = currentCursor.Position;
    pointerSize = str2double(handles.pointerSize.String);
       
    colorOrder = get(gca,'ColorOrder');
    currentColor = colorOrder(rem(length(findall(gca,'Type','rectangle')),size(colorOrder,1))+1,:);
    
    Frame = [cPos(1)-pointerSize cPos(2)-pointerSize pointerSize*2 pointerSize*2];
    hold(handles.FrameImage,'on');
    rectangle('Position',Frame,'linewidth',2,'edgecolor',currentColor,'parent',handles.FrameImage);
    hold(handles.FrameImage,'off');
    datacursormode off
    hObject.Value = false;
        
    % get data from frame stack
    data = handles.UserData.Frames;
    rotAngle = str2double(handles.rotAngel.String); %angle of rotation in degrees.
    if ~isinf(rotAngle) % Rotate image
        data = imrotate(data,rotAngle,'crop');
        mask = ~imrotate(true(size(data)),rotAngle,'crop');
        data(mask) = NaN;
    end
    dSize = size(data);
    mask = true(dSize(1:2));
    [yInd,xInd] = meshgrid(cPos(1)-pointerSize:cPos(1)+pointerSize,cPos(2)-pointerSize:cPos(2)+pointerSize); %isolate index for selected area
    
    %make sure index is not out of the image area
    ind = yInd>dSize(2) | xInd>dSize(1) | yInd<=0 | xInd<=0;
    yInd(ind) = [];xInd(ind) = []; 
    mask(xInd,yInd) = false;
    
    data = reshape(data,[numel(mask),dSize(length(size(mask))+1:end)]); %merge x and y dimension based on mask size and remaining dimensions.
    mask = reshape(mask,numel(mask),1); %reshape mask to vector
    data(mask,:) = []; %delete selected pixels
    
    dSize = size(data);
    amean = nanmean(data);
    astd = nanstd(double(data))/sqrt(dSize(1));
    if handles.showError.Value
        fill([1:dSize(2) dSize(2):-1:1],[amean+astd fliplr(amean-astd)], currentColor, 'FaceAlpha', 0.5,'linestyle','none', 'parent', handles.FrameBuffer);
    else
        fill([1:dSize(2) dSize(2):-1:1],[amean+astd fliplr(amean-astd)], currentColor, 'FaceAlpha', 0,'linestyle','none', 'parent', handles.FrameBuffer);
    end
    hold(handles.FrameBuffer,'on');
    plot(amean,'linewidth',1,'color',currentColor, 'parent', handles.FrameBuffer);
    set(handles.FrameBuffer,'xlim',[0 length(amean)]);
    axis(handles.FrameBuffer,'square');
    
%     handles.FrameBuffer.XTick = [];
%     handles.FrameBuffer.YTick = [];
    
    else
        disp('Invalid pointer target. Click on the image to plot corresponding traces.')
    end
    pause(0.5);
    delete(findall(handles.FrameImage,'Type','hggroup','HandleVisibility','off'));
    delete(findall(handles.FrameBuffer,'Type','hggroup','HandleVisibility','off'));
else
    disp('Traces can only be plotted if more than one frame is provided.')
end

% --- Executes on button press in CurrentFrameBuffer.
function CurrentFrameBuffer_Callback(hObject, eventdata, handles)
% hObject    handle to CurrentFrameBuffer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in clearTrace.
function clearTrace_Callback(hObject, eventdata, handles)
% hObject    handle to clearTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


delete(findall(handles.FrameImage.Children,'Type','rectangle'));
cla(handles.FrameBuffer);


% --- Executes on mouse press over figure background.
function ShowStack_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to ShowStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close ShowStack.
function ShowStack_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to ShowStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if handles.UserData.Close
    delete(handles.ShowStack)
else
    handles.UserData.Close = true; %leave figure but give a flag to output function to close after output has been provided
    guidata(hObject, handles);
    uiresume(hObject);
end


% --- Executes on selection change in pickColormap.
function pickColormap_Callback(hObject, eventdata, handles)
% hObject    handle to pickColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pickColormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pickColormap

cMap = str2func(strtrim(handles.pickColormap.String{handles.pickColormap.Value}));
C = cMap(256);
colormap(handles.FrameImage,C);


% --- Executes during object creation, after setting all properties.
function pickColormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pickColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ScaleColor.
function ScaleColor_Callback(hObject, eventdata, handles)
% hObject    handle to ScaleColor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if hObject.Value
    fLength = round(str2double(get(handles.FrameSmth,'string'))); %determines degree of smoothing.
    data = ApplyFilter2(squeeze(handles.UserData.Frames(:,:,get(handles.FrameSlider,'Value'))),fLength,str2double(handles.FrameSigma.String),get(handles.FilterType,'value')); %smooth current frame
    
    caxis(handles.FrameImage,[min(data(:)) max(data(:))]);
    handles.MinColor.String = num2str(round(double(min(data(:))),4));
    handles.MaxColor.String = num2str(round(double(max(data(:))),4));
end

% --- Executes on button press in SaveMovie.
function SaveMovie_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMovie (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fName,fPath] = uiputfile('*.avi', 'Enter filename');

cMap = str2func(strtrim(handles.pickColormap.String{handles.pickColormap.Value}));
C = cMap(256);
Range = [str2double(handles.MinColor.String) str2double(handles.MaxColor.String)];

% v = VideoWriter([fPath fName],'indexed avi');
% v.Colormap = C;
v = VideoWriter([fPath fName]);
v.FrameRate = str2double(inputdlg({'Enter the framerate for movie file'},'',1));
open(v);

for iFrames = 1:size(handles.UserData.Frames,3)
    
    pic = ApplyFilter2(squeeze(handles.UserData.Frames(:,:,iFrames)), ... %smooth current frame
        round(str2double(get(handles.FrameSmth,'string'))), ...
        str2double(handles.FrameSigma.String),get(handles.FilterType,'value'));
    
    rotAngle = str2double(handles.rotAngel.String); %angle of rotation in degrees.
    if ~isinf(rotAngle) % Rotate image
        pic = imrotate(pic,rotAngle,'crop');
    end
    pic(isnan(pic)) = 0;
    
    if handles.UseThreshold.Value
        thresh = str2double(handles.ImgThresh.String);
        if thresh >= 0
            ind = uint8(imfill(pic > nanstd(pic(:))*thresh, 'holes')); %threshold input matrix and make sure there are no holes in separate areas
        elseif thresh < 0
            ind = uint8(imfill(pic < nanstd(pic(:))*thresh, 'holes')); %threshold input matrix and make sure there are no holes in separate areas
        end
        pic = bsxfun(@times, double(pic), double(ind)); %set values below thresh to 0. Will come out as black.
    end

    pic = mat2gray(pic,Range);
    pic = gray2ind(pic,256);
    pic = ind2rgb(pic, C);
    writeVideo(v,pic)
    
end
close(v)    


% --- Executes on button press in UseThreshold.
function UseThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to UseThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

drawImage(handles,handles.FrameImage,handles.FrameSlider);


function ImgThresh_Callback(hObject, eventdata, handles)
% hObject    handle to ImgThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImgThresh as text
%        str2double(get(hObject,'String')) returns contents of ImgThresh as a double

cVal = str2double(hObject.String);
if isnan(cVal)
    hObject.String = 0;
    disp([get(hObject,'string') ' is not a valid input for threshold. Set back to 0.'])
else
    UseThreshold_Callback([], [], handles); %re-draw images.
end

% --- Executes during object creation, after setting all properties.
function ImgThresh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImgThresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in releaseFigure.
function releaseFigure_Callback(hObject, eventdata, handles)
% hObject    handle to releaseFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume(gcbf); %return control to matlab

% --- Executes on button press in selectArea.
function selectArea_Callback(hObject, eventdata, handles)
% hObject    handle to selectArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

threshold = str2double(handles.ImgThresh.String); %threshold in standard deviation units. This needs to be at least

if threshold ~= 0
    datacursormode on
    cursorMode = datacursormode(handles.ShowStack);
    currentCursor = getCursorInfo(cursorMode);
    while isempty(currentCursor)
        drawnow;
        currentCursor = getCursorInfo(cursorMode);
        pause(0.01);
    end
    
    cPos = currentCursor.Position;
    colorOrder = repmat(reshape([0:1/4:1;1:-1/4:0],10,1),1,3); %color order is going from white to black
    currentColor = colorOrder(rem(length(findall(gca,'Type','line')),size(colorOrder,1))+1,:);
    axLabel = currentCursor.Target.Parent.Tag; %get label of selected axis
    
    % get frame data, based on which axis was selected
    fLength = round(str2double(get(handles.FrameSmth,'string'))); %determines degree of smoothing.
    data = ApplyFilter2(squeeze(handles.UserData.Frames(:,:,get(handles.FrameSlider,'Value'))),fLength,str2double(handles.FrameSigma.String),get(handles.FilterType,'value')); %smooth current frame
    
    rotAngle = str2double(handles.rotAngel.String); %angle of rotation in degrees.
    if ~isinf(rotAngle) % Rotate image
        data = imrotate(data,rotAngle,'crop');
        mask = ~imrotate(true(size(data)),rotAngle,'crop');
        data(mask) = NaN;
    end
    
    if threshold >= 0
        ind = imfill(data > nanstd(data(:))*threshold, 'holes'); %threshold input matrix and make sure there are no holes in separate areas
    elseif threshold < 0
        ind = imfill(data < nanstd(data(:))*threshold, 'holes'); %threshold input matrix and make sure there are no holes in separate areas
    end
    
    if ~isempty(ind)
        foundAreas = bwlabel(ind, 8); %find independent areas
        cArea = foundAreas(cPos(2),cPos(1)); %get selected area
        areaOutline = bwboundaries(foundAreas == cArea); %outline of selected area
        
        % plot to both axis
        hold(handles.FrameImage,'on');
        plot(handles.FrameImage,areaOutline{1}(:,2),areaOutline{1}(:,1),'linewidth',4,'color',currentColor)
        hold(handles.FrameImage,'off');
    end
    datacursormode off
        
    pause(0.5);
    delete(findall(handles.FrameImage,'Type','hggroup','HandleVisibility','off'));
end


% --- Executes on selection change in areaProp.
function areaProp_Callback(hObject, eventdata, handles)
% hObject    handle to areaProp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns areaProp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from areaProp


% --- Executes during object creation, after setting all properties.
function areaProp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to areaProp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in clearAreas.
function clearAreas_Callback(hObject, eventdata, handles)
% hObject    handle to clearAreas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(findall(handles.FrameImage.Children,'Type','line')); %delete all lines from images - this should take care of area outlines


% --- Executes on button press in closeFigure.
function closeFigure_Callback(hObject, eventdata, handles)
% hObject    handle to closeFigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ShowStack_CloseRequestFcn(handles.ShowStack, [], handles); %close figure


function rotAngel_Callback(hObject, eventdata, handles)
% hObject    handle to rotAngel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of rotAngel as text
%        str2double(get(hObject,'String')) returns contents of rotAngel as a double

if isnan(str2double(get(hObject,'string')))
    set(hObject,'string','0');
    disp([get(hObject,'string') ' is not a valid input for rotation angle. Set back to 0.'])
else
    FilterType_Callback([],[],handles);
end

% --- Executes during object creation, after setting all properties.
function rotAngel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to rotAngel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%% Additional functions
function DataOut = ApplyFilter2(DataIn,fLength,sigma,type)
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
DataOut = padarray(double(DataOut),[1 1],NaN,'both'); %add NaNs to image to see its edge against zero-padded background when rotating.


function drawImage(handles,currAxis,cSlider)
% Short function to draw image to selected axis. 
% handles are figure handles. currAxis is the current axis (frameImage or
% frame buffer) and cSlider is the according frame slider.

fLength = str2double(get(handles.FrameSmth,'String'));
data = ApplyFilter2(squeeze(handles.UserData.Frames(:,:,get(cSlider,'Value'))),fLength,str2double(handles.FrameSigma.String),get(handles.FilterType,'value')); %smooth current frame

hold(currAxis,'on')
delete(findall(currAxis.Children,'Type','image'));
rotAngle = str2double(handles.rotAngel.String); %angle of rotation in degrees. 
if ~isinf(rotAngle) % Rotate image
    data = imrotate(data,rotAngle,'crop');
    mask = ~imrotate(true(size(data)),rotAngle,'crop');
    data(mask) = NaN;
end
cImg = imshow(data,[str2double(get(handles.MinColor,'String')) str2double(get(handles.MaxColor,'String'))],'parent',currAxis); %create image object for preview
colormap(handles.FrameImage,strtrim(handles.pickColormap.String{handles.pickColormap.Value}));
set(cImg,'AlphaData',~isnan(data)); %make NaNs transparent.
uistack(cImg,'bottom'); %make this bottom of the stack so squares remain visible
hold(currAxis,'off')

if handles.ScaleColor.Value %adjust color range if set to auto scale
    if strcmp(currAxis.Tag,'FrameImage')
        cRange = [min(data(:)) max(data(:))];
        if ~any(isnan(cRange))
            caxis(currAxis,cRange);
            handles.MinColor.String = num2str(round(double(min(data(:))),4));
            handles.MaxColor.String = num2str(round(double(max(data(:))),4));
        end
    end
end

if handles.UseThreshold.Value
    ind = thresholdImage(data,str2double(handles.ImgThresh.String));
    set(cImg,'AlphaData',ind); %make values below transparent. Since there is no image below, pixels should appear white.
end


function ind = thresholdImage(data,thresh)
% quick function to do the thresholding

if thresh >= 0
    ind = imfill(data > nanstd(data(:))*thresh, 'holes'); %threshold input matrix and make sure there are no holes in separate areas
elseif thresh < 0
    ind = imfill(data < nanstd(data(:))*thresh, 'holes'); %threshold input matrix and make sure there are no holes in separate areas
end


%additional colormap
function map = colormap_blueblackred(m)
if nargin < 1
   m = size(get(gcf,'colormap'),1);
end

map=[1 1 0; 1 0.96 0; 1 0.92 0; 1 0.88 0; 1 0.84 0; 1 0.80 0; 1 0.76 0; 1 0.72 0; 1 0.68 0; 1 0.64 0; 1 0.60 0; 1 0.56 0; 1 0.52 0; 1 0.48 0; 1 0.44 0; 1 0.40 0;  ...
    1 0.36 0; 1 0.32 0; 1 0.28 0; 1 0.24 0; 1 0.20 0; 1 0.16 0; 1 0.12 0; 1 0.08 0; 1 0.04 0;
    1 0 0; 0.96 0 0; 0.92 0 0; 0.88 0 0; 0.84 0 0; 0.80 0 0; 0.76 0 0; 0.72 0 0; 0.68 0 0; 0.64 0 0; 0.60 0 0; 0.56 0 0; 0.52 0 0; 0.48 0 0; 0.44 0 0; 0.40 0 0;  ...
    0.36 0 0; 0.32 0 0; 0.28 0 0; 0.24 0 0; 0.20 0 0; 0.16 0 0; 0.12 0 0; 0.08 0 0; 0.04 0 0; 0 0 0;                                   ...
    0 0 0.04;  0 0 0.08; 0 0 0.12; 0 0 0.16; 0 0 0.20; 0 0 0.24; 0 0 0.28; 0 0 0.32; 0 0 0.36; 0 0 0.40; 0 0 0.44; 0 0 0.48; 0 0 0.52; ...
    0 0 0.56; 0 0 0.60; 0 0 0.64; 0 0 0.68; 0 0 0.72; 0 0 0.76; 0 0 0.80; 0 0 0.84; 0 0 0.88; 0 0 0.92; 0 0 0.96; 0 0 1; ...
    0 0.04 1;  0 0.08 1; 0 0.12 1; 0 0.16 1; 0 0.20 1; 0 0.24 1; 0 0.28 1; 0 0.32 1; 0 0.36 1; 0 0.40 1; 0 0.44 1; 0 0.48 1; 0 0.52 1; ...
    0 0.56 1; 0 0.60 1; 0 0.64 1; 0 0.68 1; 0 0.72 1; 0 0.76 1; 0 0.80 1; 0 0.84 1; 0 0.88 1; 0 0.92 1; 0 0.96 1; 0 1 1];

map = imresize(map,[m,3]);
map = map - min(map(:));map = map ./ max(map(:));
map = map(end:-1:1,:);


% --- Executes on button press in showError.
function showError_Callback(hObject, eventdata, handles)
% hObject    handle to showError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showError

patches = findall(handles.FrameBuffer.Children,'Type','Patch');
for iPatch = 1 : size(patches,1)
    if handles.showError.Value
        patches(iPatch).FaceAlpha = 0.5;
    else
        patches(iPatch).FaceAlpha = 0;
    end
end

% --- Executes on button press in showOutline.
function showOutline_Callback(hObject, eventdata, handles)
% hObject    handle to showOutline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showOutline




% --- Executes on button press in loadFrame.
function loadFrame_Callback(hObject, eventdata, handles)
% hObject    handle to loadFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fName, fPath] = uigetfile(pwd,'Select single frame');

data = load([fPath fName]);
cField = fieldnames(data);

if length(size(data.(cField{1}))) == 2
    hold(handles.FrameBuffer,'off');
    delete(findall(handles.FrameImage.Children,'Type','rectangle'));
    cla(handles.FrameBuffer);
    imagesc(data.(cField{1}),'parent',handles.FrameBuffer);
    axis(handles.FrameBuffer,'image');
    colormap(handles.FrameBuffer,'gray');
    view(0,90);
    handles.FrameBuffer.XTick = [];
    handles.FrameBuffer.YTick = [];
else
    disp('Selected file must contain a single 2D matrix to plot');
end


% --- Executes on button press in loadOutline.
function loadOutline_Callback(hObject, eventdata, handles)
% hObject    handle to loadOutline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
