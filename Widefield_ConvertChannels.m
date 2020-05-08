opts.fPath = 'C:\data\WidefieldImager\Animals\mSM34\PhaseMap\22-Mar-2017';
opts.fName = 'Frames';
opts.plotChans = true;
opts.trigLine = NaN;

for iTrials = 6:60
    
%     [blueData,blueTimes,hemoData,hemoTimes] = Widefield_CheckChannels(opts,1);
    [blueData,blueTimes] = Widefield_CheckChannels(opts,iTrials);
    
    cFile = [opts.fPath '\blueFrames_' num2str(iTrials) '.dat'];
    Widefield_SaveData(cFile,blueData,blueTimes); %write new file for aligned blue channel.
    
%     cFile = [opts.fPath '\hemoFrames_' num2str(iTrials) '.dat'];
%     Widefield_SaveData(cFile,hemoData,hemoTimes); %write new file for aligned hemo channel.
    
    cFile = [opts.fPath '\' opts.fName '_' num2str(iTrials) '.dat']; %current file to be read
    delete(cFile);
   
end