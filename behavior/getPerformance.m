clear,clc,close all

%% load data object (one session)

datapth = '/Users/Munib/Documents/Economo-Lab/code/data'; % data path
datafn = 'data_structure_JEB7_2021-04-29'; % file name

temp = load(fullfile(datapth,datafn));
obj = temp.obj;
clear temp;

%% performance

% RIGHT trials where animal did not lick early and was not provided free water
try
    % for new data objects, autowater is a logical array
    nTrials = sum(obj.bp.R & ~obj.bp.early & ~obj.bp.autowater);
    nTrialsCorrect = sum(obj.bp.R & obj.bp.hit & ~obj.bp.early & ~obj.bp.autowater);
catch
    % for old data objects, autowater is a struct
    nTrials = sum(obj.bp.R & ~obj.bp.early & obj.bp.autowater.nums==2);
    nTrialsCorrect = sum(obj.bp.R & obj.bp.hit & ~obj.bp.early & obj.bp.autowater.nums==2);
end

rightPerformance = nTrialsCorrect / nTrials * 100; % STORE THIS VALUE IN DATABASE

% LEFT trials where animal did not lick early and was not provided free water
try
    % for new data objects, autowater is a logical array
    nTrials = sum(obj.bp.L & ~obj.bp.early & ~obj.bp.autowater);
    nTrialsCorrect = sum(obj.bp.L & obj.bp.hit & ~obj.bp.early & ~obj.bp.autowater);
catch
    % for old data objects, autowater is a struct
    nTrials = sum(obj.bp.L & ~obj.bp.early & obj.bp.autowater.nums==2);
    nTrialsCorrect = sum(obj.bp.L & obj.bp.hit & ~obj.bp.early & obj.bp.autowater.nums==2);
end

leftPerformance = nTrialsCorrect / nTrials * 100; % ALSO STORE THIS VALUE IN DATABASE










