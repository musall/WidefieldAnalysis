function Widefield_RawToMpeg(fPath,skipCheck)
% code to convert widefield raw data from binary (.dat) to .mj2 compressed video files

if ~exist('skipCheck','var')
    skipCheck = false; %only convert data if Vc, blueV and hemoV are present
end

if fPath(end) ~= filesep
    fPath = [fPath filesep];
end

%% folders that contain animal video data
animals = dir(fPath);
animals = animals([animals.isdir] & ~strncmpi('.', {animals.name}, 1));

%% loop through animal data
for iAnimals = 1:length({animals.name})
    
    aPath = [animals(iAnimals).name filesep];
    paradigms = dir([fPath aPath]);
    paradigms = paradigms([paradigms.isdir] & ~strncmpi('.', {paradigms.name}, 1));
    
    %% loop through paradigms
    for iParadigms = 1:length({paradigms.name})
        
        pPath = [paradigms(iParadigms).name filesep];
        sessions = dir([fPath aPath pPath]);
        sessions = sessions([sessions.isdir] & ~strncmpi('.', {sessions.name}, 1));
        
        %% loop through sessions
        for iSessions = 1:length({sessions.name})
            
            sPath = [sessions(iSessions).name filesep];
            files = dir([fPath aPath pPath sPath  'Frames_*.dat']);
            
            VcCheck = exist([fPath aPath pPath sPath 'Vc.mat'],'file') == 2;
            blueVCheck = exist([fPath aPath pPath sPath 'blueV.mat'],'file') == 2;
            hemoVCheck = exist([fPath aPath pPath sPath 'hemoV.mat'],'file') == 2;
            
            if ~isempty(files)
                if skipCheck || (VcCheck && blueVCheck && hemoVCheck)
                    disp(['Found Vc in experiment: ' animals(iAnimals).name ' - ' sPath(1:end-1) '. Converting raw data to .mp4.']);
                    tic
                    for iFiles = 1:length(files)
                        
                        cFile = [fPath aPath pPath sPath files(iFiles).name]; %current file to be read
                        [frameTimes,Data] = Widefield_LoadData(cFile,'Frames'); %load video data
                        frameTimes = frameTimes(1:end-length(size(Data))); %extract frame times from header and convert to millisecond timestamps
                        
                        cFile = strrep(cFile,'.dat','.mat');
                        save(strrep(cFile,'Frames','frameTimes'),'frameTimes'); %save frametimes
                        
                        if max(Data(:)) > 255 %data is not int8
                            if max(Data(:)) > 65535 %data is not int16
                                error('Unknown data type');
                            else
                                Data = mat2gray(Data,[0 65535]);
                            end
                        else
                            Data = mat2gray(Data,[0 255]);
                        end
                        
                        Data = padarray(Data, 8 - [rem(size(Data,1),8), rem(size(Data,2),8)], 0, 'post'); %pad array to be dividable by 8
                        Data = reshape(Data, size(Data,1),size(Data,2),1,[]);
                        
                        v = VideoWriter(strrep(cFile, '.mat', ''), 'MPEG-4'); %save as compressed video file
                        v.Quality = 100;
                        open(v);
                        if size(Data,4) > 500
                            for iFrames = 1 : size(Data,4)
                                writeVideo(v,Data(:,:,:,iFrames));
                            end
                        else
                            writeVideo(v,Data);
                        end
                        close(v);

                        cFile = strrep(cFile,'mat','dat');
                        delete(cFile);
                        
                    end
                    toc
                end
            end
        end
    end
end