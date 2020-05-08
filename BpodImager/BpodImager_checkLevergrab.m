function [grabOn,grabRel,tapOn,tapRel] = BpodImager_checkLevergrab(tapDur,postStimDur,grabs,release)

grabOn = [];
grabRel = [];
tapOn = [];
tapRel = [];

if ~isempty(grabs)
    Cnt = 0;
    grabCnt = 0;
    tapCnt = 0;
    
    while true
        Cnt = Cnt +1;
        if length(grabs) < Cnt %no more grabs
            break;
        end
        
        idx = find(release > grabs(Cnt),1);
        if ~isempty(idx) %if there is a release that follows current grab
            cGrab = release(idx) - grabs(Cnt); %duration of current grab
        else
            cGrab = postStimDur - grabs(Cnt);
        end
        
        %% more grabs available - check start time of next grab and merge if they are very close in time
        if length(grabs) > Cnt && ~isempty(idx)
            while (grabs(Cnt+1) - release(idx)) <= 0.5 %time diff between grabs is less than 50 ms. Merge with next grab.
                
                release(idx) = []; %delete current release
                grabs(Cnt+1) = []; %delete next grab
                
                idx = find(release > grabs(Cnt),1); %next release
                if ~isempty(idx) %if there is a release that follows current grab
                    cGrab = release(idx) - grabs(Cnt); %duration of current grab
                else
                    cGrab = postStimDur - grabs(Cnt);
                end
                                
                if isempty(idx) || length(grabs) <= Cnt %no more grabs/releases
                    break;
                end

            end
        end
        
        %% check if current grab is grab or tap
        if cGrab <= tapDur
            tapCnt = tapCnt + 1;
            tapOn(tapCnt) = grabs(Cnt);
            if ~isempty(idx)
                tapRel(tapCnt) = idx;
            end
        else
            grabCnt = grabCnt + 1;
            grabOn(grabCnt) = grabs(Cnt);
            if ~isempty(idx)
                grabRel(grabCnt) = idx;
            end
        end
    end
end



