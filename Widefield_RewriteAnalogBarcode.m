function Widefield_RewriteAnalogBarcode(opts,firstTrial,lastTrial,fileIdx)
% code to replace trialces in the analog barcode line. This is to repair
% session with multiple behavior files or if barcodes were not recorded correctly.
% firstTrial and lastTrial denote a range of trials from which analog files
% will be loaded and barcode will be writen to match the trialNr as defined
% by their fileName. This is to allow combining analog files when behavior
% had to be restarted.
%
% If fileIdx is provided, the code expects a 2xtrials matrix where the
% first row indicates the trialNumber of the analog file and the second row
% the fileNr that should be written.

if ~exist('fileIdx','var') || isempty(fileIdx)
    fileIdx = repmat(firstTrial : lastTrial,2,1);
end

for iTrials = 1 : size(fileIdx,2)
    
    cFile = [opts.fPath filesep 'Analog_' num2str(fileIdx(1,iTrials)) '.dat']; %analog file to be read
    [header,Analog] = Widefield_LoadData(cFile,'Analog'); %load analog data
    
    temp = ones(1,size(Analog,2)) * 3500;
    temp(find(diff(double(Analog(opts.barcodeLine,:)) < 1000),1):end) = 0;
    
    code = encode2of5(fileIdx(2,iTrials));
    codeModuleDurs = [3 6]; %Durations for each module of the trial code sent over the TTL line
    
    Cnt = 1;
    barSeq = zeros(1,sum(code(:)) * 2);
    for iCode = 1:size(code,2)
        
        barSeq(Cnt : Cnt + codeModuleDurs(code(1,iCode)) -1) = 3500;
        barSeq(Cnt + codeModuleDurs(code(1,iCode)) : Cnt + codeModuleDurs(code(1,iCode)) + codeModuleDurs(code(2,iCode)) -1) = 0;
        
        Cnt = Cnt + sum(codeModuleDurs(code(:,iCode)));
        
    end
    
    ind = find(diff(double(Analog(opts.barcodeLine,:)) < 1000),1) + 50;
    if isempty(ind)
        Analog(opts.barcodeLine,:) = zeros(1,size(Analog,2));
        Analog(opts.barcodeLine,1:200) = ones(1,200)*3500;
        ind = 250;
    end
    temp(ind + 1: ind + length(barSeq)) = barSeq;
    Analog(opts.barcodeLine,:) = uint16(temp);
    
    fid = fopen(cFile, 'wb'); %open binary stimulus file
    fwrite(fid,3,'double'); %indicate number of single values in the header
    fwrite(fid,header,'double'); %write number of recorded analog channels + timestamps
    fwrite(fid,Analog,'uint16');
    fclose(fid);
    
end
