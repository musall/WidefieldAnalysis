%% combine behavior
firstPath = 'U:\space_managed_data\BehaviorVideo\Fez10\SpatialDisc\Session Data\Fez10_SpatialDisc_Mar13_2019_Session1';
secondPath = 'U:\space_managed_data\BehaviorVideo\Fez10\SpatialDisc\Session Data\Fez10_SpatialDisc_Mar13_2019_Session2';
idx1 = 1:223;
idx2 = 1:189;
gapTrials = 1; %nr of imaging trials that were recorded between behavioral sessions

[a, b] = fileparts(firstPath);
load([a filesep b '.mat'])
bhv1 = selectBehaviorTrials(SessionData,idx1);
if sum(SessionData.DidNotLever) < gapTrials
    error('not enough unperformed trials that can be used as filler');
end
gapBhv = selectBehaviorTrials(SessionData,find(SessionData.DidNotLever,gapTrials)); %add empty trials (this is to compensate trials where bpod broke)
bhv1 = appendBehavior(bhv1,gapBhv); %combine bhv data


[a, b] = fileparts(secondPath);
load([a filesep b '.mat'])
bhv2 = selectBehaviorTrials(SessionData,idx2);

SessionData = appendBehavior(bhv1,bhv2); %combine bhv data
SessionData.nTrials = idx1(end)+idx2(end)+gapTrials;

temp = findstr(b,'_');
cFile = [a filesep b(1:temp(end)) 'Combined'];
save(cFile, 'SessionData');

%% rename behavior video
Widefield_RenameBehaviorVideo(firstPath,secondPath,idx1(end)+gapTrials)

%% rewrite analog data
opts.fPath = 'U:\smusall\BpodImager\Animals\Fez10\SpatialDisc\13-Mar-2019';
opts.barcodeLine = 7; %line for Bpod trialID
Widefield_RewriteAnalogBarcode(opts,377,382)
