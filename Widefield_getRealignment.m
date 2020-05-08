function newVc = Widefield_getRealignment(Vc, cBhv, segFrames, opts)
% align imaging data using Session data
% this will align to time of handle grab and stimulus onset

newVc = NaN(size(Vc,1), segFrames(end), size(Vc,3), 'single'); %new Vc to capture max duration of each segment
for iTrials = 1 : size(Vc,3)

    % get indices for current trial
    stimOn = cBhv.RawEvents.Trial{iTrials}.Events.Wire3High; %time of stimulus onset - measured from soundcard
    handleOn = [reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal1',1,[]) ...
        reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal2',1,[]) ...
        reshape(cBhv.RawEvents.Trial{iTrials}.States.WaitForAnimal3',1,[])];
    
    clear cIdx
    cIdx(1) = handleOn(find(handleOn == cBhv.RawEvents.Trial{iTrials}.States.WaitForCam(1))-1); %find start of lever state that triggered stimulus onset
    cIdx(2) = stimOn;
    cIdx = floor((cIdx - stimOn + opts.preStim) * opts.frameRate); %convert to frames. This is the last frame of each segment.
    cIdx(end + 1) = size(Vc,2);
    
    if segFrames(1) >= cIdx(1)
        newVc(:, segFrames(1) - cIdx(1) + 1 : segFrames(1), iTrials) = Vc(:, 1 : cIdx(1), iTrials); % baseline
    elseif segFrames(1) < cIdx(1)
        newVc(:, 1 : segFrames(1), iTrials) = Vc(:, cIdx(1) - segFrames(1) + 1 : cIdx(1), iTrials); % baseline
    end
        
    newVc(:, segFrames(1) + 1 : segFrames(1) + (diff(cIdx(1:2))), iTrials) = Vc(:, cIdx(1) + 1 : cIdx(2), iTrials); %handle period
    newVc(:, segFrames(2) + 1 : end, iTrials) = Vc(:, cIdx(2) + 1 : cIdx(2) + diff(segFrames(2:3)), iTrials); %stimulus period
    
end