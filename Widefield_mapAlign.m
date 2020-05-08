function data = Widefield_mapAlign(data,opts)
% short code to align widefield data. 
% Straightens image and shifts so bregma is at the center point. 

dSize = size(data);
if length(dSize) > 2
    temp = ~imrotate(true(dSize(1:2)), opts.rotAngle,'crop');
    data = imrotate(data, opts.rotAngle,'crop');
    data = reshape(data,[],dSize(3));
    if ~islogical(data)
        data(temp(:),:) = NaN;
    else
        data(temp(:),:) = false; 
    end
    data = reshape(data,dSize);
else
    data = imrotate(data, opts.rotAngle,'crop');
    if ~islogical(data)
        data(~imrotate(true(dSize), opts.rotAngle,'crop')) = NaN;
    else
        data(~imrotate(true(dSize), opts.rotAngle,'crop')) = false;
    end
end

%shift matrices to put bregma in the center
shiftY = round((dSize(1) - (opts.bregma(2)*2)) / 2);
shiftX = round((dSize(2) - (opts.bregma(1)*2)) / 2);
if ~islogical(data)
    data = arrayShift(data,[shiftX shiftY],NaN);
else
    data = logical(arrayShift(data,[shiftX shiftY],false));
end

% cut to size if NaN mask was used
if sum(~isnan(data(:))) ~= numel(data) && sum(~isnan(data(:))) > 0
    mask = ~isnan(data(:,:,1));
    xCut = find(sum(mask,1) > 0);
    xCut = xCut([1 end]);
    xCut = min([xCut(1) size(mask,2)-xCut(end)]);
    
    yCut = find(sum(mask,2) > 0);
    yCut = yCut([1 end]);
    yCut = min([yCut(1) size(mask,1)-yCut(end)]);
    
    data = data(yCut:end-yCut+1,xCut:end-xCut+1,:);
    dSize(1) = size(data,1);dSize(2) = size(data,2);
    data = reshape(data,dSize);
end
