clear,clc,close all

% add folders in current directory to the path
% this way you have access to functions in funcs and utils directories
addpath(genpath(pwd))


%% SET PARAMS
params.alignEvent          = 'jawOnset'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};          % right hits, no stim, aw on
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};          % left hits, no stim, aw on
% params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};        % error right, no stim, aw off
% params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};        % error left, no stim, aw off
% params.condition(end+1) = {'~hit&~miss&~stim.enable&~autowater&~early'};    % ignore
% params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};           % hit 2afc
% params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};            % hit aw

% specify probe number and areas to load and process data
params.probe(1) = 1;
params.probeArea{1} = 'xxx';
% 
% params.probe(end+1) = 2;
% params.probeArea{end+1} = 'ALM';

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

%% SET METADATA
% experiment meta data
meta.datapth = '/Users/Munib/Documents/Economo-Lab/code/data';
% meta.datapth = '/Volumes/MUNIB_SSD/Experiments';
meta.anm = 'JGR2'; % 'JEB7'  'EKH3'  'JGR2'
meta.date = '2021-11-16'; % '2021-04-29'  '2021-08-11'  '2021-11-16'
meta.datafn = findDataFn(meta);

%% LOAD AND PROCESS DATA
dat = load(fullfile(meta.datapth,meta.datafn));
obj = dat.obj;

for prbnum = 1:numel(params.probe)
    disp('______________________________________________________')
    disp(['Processing data for probe ' num2str(prbnum)])
    disp(' ')
    [params,obj] = processData(obj,params,prbnum);
end

% if only one probe, clean up so all previous code works
if numel(params.probe)==1
    [obj,params] = oneProbeTrim(obj,params);
end

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')

%% find early move trials, movement and non-movement times using motion energy

% TODO - make sure go cue is aligned (might have to add 0.5 sec to time
% vector from me.data
% - also return ix of movement times based on obj.time

params.moveThresh = 0.15; % fraction of time in delay period that if labeled as moving, labels trial as a early move trial
[meta,mov] = assignEarlyTrials(obj,meta,params);



















