clear,clc,close all

% add folders in current directory to the path
% this way you have access to functions in funcs and utils directories
addpath(genpath(pwd))


%% SET METADATA (ANIMALS AND RECORDING SESSIONS) AND PROBE INFO

% experiment meta data (assumes all data objs in same directory)
meta.datapth = '/Users/Munib/Documents/Economo-Lab/code/data';
% meta.datapth = '/Volumes/MUNIB_SSD/Experiments';

meta.anm = {'EKH3';...
            'JEB7'};

meta.date = { {'2021-08-11'};...
              {'2021-04-29'} };
          
          
probeInfo.probe = {[1 2];...
                   1};
probeInfo.probeArea = { {'IRN','ALM'};...
                        {'ALM'} };

%% LOAD AND PROCESS DATA

% NEED TO COME UP WITH A WAY OF SEPARATING META ABOVE FROM WHAT I NEED FOR
% EACH ANM/SESSION


% find data file names for each animal and recording date
% then load and process
max_size = findMaxDates(meta.date);
meta.datafn = cell(numel(meta.anm),max_size);
ct = 1;
for i = 1:numel(meta.anm)
    for j = 1:numel(meta.date{i})
        
        meta.datafn{i,j} = findDataFnMultiSession(meta,i,j);
        [params{i,j} obj{i,j}] = main_multiSession(probeInfo,meta)
        ct = ct + 1;
    end
end


%%

function max_size = findMaxDates(dates)
for i = 1:numel(dates)
    numDates(i) = numel(dates{i});
end
max_size = max(numDates);
end