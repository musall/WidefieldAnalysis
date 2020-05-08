function Widefield_CombineSessions(path1,path2)
%short function to combine two imaging sessions. This is used if bpod was
%interrupted and restarted and two sessions are created by accident.

recs = dir([path1 filesep 'Analog*']);
for iRecs = 1:length(recs)
    a = textscan(recs(iRecs).name,'%s%f%s','Delimiter','_');
    trials(iRecs) = a{2};
end
maxTrial = max(trials); clear trials
recs = dir([path2 filesep '*.dat']);

% move files
for x = 1:length(recs)
    y = textscan(recs(x).name,'%s%f%s','Delimiter','_');
    y = y{2};
    cFile = strrep(recs(x).name,num2str(y),num2str(y+maxTrial));
    movefile([path2 filesep recs(x).name],[path1 filesep cFile]);
end

%combine behavioral files. Only works if only one of them is present in each folder.
[~,thisDate] = fileparts(path1);
thisDate = textscan(thisDate,'%s%s','Delimiter','_');
thisDate = thisDate{1};
thisDate = char(datetime(thisDate,'Format','MMMdd_yyyy'));
bhvFile1 = dir([path1 filesep '*' thisDate '*.mat']);
bhvFile1 = [path1 filesep bhvFile1.name];
load(bhvFile1);
SessionData1 = SessionData;

[~,thisDate] = fileparts(path2);
thisDate = textscan(thisDate,'%s%s','Delimiter','_');
thisDate = thisDate{1};
thisDate = char(datetime(thisDate,'Format','MMMdd_yyyy'));
bhvFile2 = dir([path2 filesep '*' thisDate '*.mat']);
bhvFile2 = [path2 filesep bhvFile2.name];
load(bhvFile2);
SessionData2 = SessionData;

fPath = [path1 filesep 'SplitSessionData'];
mkdir(fPath);
[~,fName1] = fileparts(bhvFile1);
[~,fName2] = fileparts(bhvFile2);
movefile(bhvFile1,[fPath filesep fName1 '.mat']);
movefile(bhvFile2,[fPath filesep fName2 '.mat']);

SessionData = appendBehavior(SessionData1,SessionData2);
save([path1 filesep fName1 '_combined.mat'],'SessionData');


function bhv = appendBehavior(bhv,data)
% Function to collect data from behavioral files into one unified array bhv.
% Usage: bhv = appendBehavior(bhv,data)

%% get fieldnames
if isempty(bhv)
    bFields = {};
else
    bFields = fieldnames(bhv);
end
dFields = fieldnames(data);
dFields(strcmpi(dFields,'Settings')) = [];

%% cycle trough fields and add current data to bhv structure. Create new bhv entries if required.
structCases = {'RawEvents' 'RawData'}; %cases at which substructs are present. Go one level depper to correctly append them together.
for iFields = 1:size(dFields,1)
    if isstruct(data.(dFields{iFields}))
        if ismember(dFields{iFields},bFields) %existing field
            if any(strcmpi((dFields{iFields}),structCases))
                tFields = fieldnames(data.(dFields{iFields}));
                for x = 1:length(tFields)
                    bhv.(dFields{iFields}).(tFields{x}) = [bhv.(dFields{iFields}).(tFields{x}) data.(dFields{iFields}).(tFields{x})];
                end
            else
                bhv.(dFields{iFields}) = [bhv.(dFields{iFields}) data.(dFields{iFields})];
            end
        else %new field in data
            bhv.(dFields{iFields}){1} = data.(dFields{iFields});
        end
    else
        if length(data.(dFields{iFields})) > 1 %vector or matrix
            if ischar(data.(dFields{iFields})) % carry strings in a cell container
                temp = {data.(dFields{iFields})};
            else
                if any(size(data.(dFields{iFields})) == 1) && length(data.(dFields{iFields})) >= sum(data.nTrials) %check if vector with at least nTrials entries
                    temp = data.(dFields{iFields})(1:sum(data.nTrials)); %keep nTrials entries
                else
                    if iscell(data.(dFields{iFields}))
                        temp = {data.(dFields{iFields})}; %cells are usually settings. keep in one cell per data file.
                    else
                        temp = data.(dFields{iFields}); %use all values and append together into larger matrix.
                    end
                end
            end
        else %single value
            temp = data.(dFields{iFields});
        end
        
        if ~isobject(temp) % don't append objects into larget structure
            if size(temp,2) == 1; temp = temp'; end %column vectors are transposed to rows
            
            if any(ismember(dFields,bFields)) %existing field
                bhv.(dFields{iFields}) = [bhv.(dFields{iFields}) temp];
            else %new field in data
                bhv.(dFields{iFields}) = temp;
            end
        end
    end 
end