clear,clc,close all

if ispc
    pth = 'C:\Code\activity_modes';
elseif ismac
    pth = '/Users/Munib/Documents/Economo-Lab/code/activity_modes';
end

addpath(genpath(pwd))


%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

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
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};        % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};        % error left, no stim, aw off
params.condition(end+1) = {'~hit&~miss&~stim.enable&~autowater&~early'};    % ignore
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};           % hit 2afc
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};            % hit aw


% set conditions used for finding activity modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %rhit, aw off 
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %lhit, aw off 
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %rmiss, aw off 
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %lmiss, aw off 
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};    % hit, aw off 

% specify probe number and areas to load and process data
params.probe(1) = 1;
params.probeArea{1} = 'ALM';

% params.probe(end+1) = 2; % for multiprobe stuff, not ready yet (TODO)
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
meta.datapth = '/Users/Munib/Documents/Economo-Lab/data';
% meta.datapth = '/Volumes/MUNIB_SSD/Experiments';
meta.anm = 'JEB7'; % 'JEB7'  'EKH3'
meta.date = '2021-04-29'; % '2021-04-29'   '2021-08-11
meta.datafn = findDataFn(meta);

%% LOAD AND PROCESS DATA
dat = load(fullfile(meta.datapth,'DataObjects',meta.anm,meta.datafn));
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

%% label move or non-move
% [obj.movix,obj.movtime] = getMoveIdx(obj,params);

%% ACTIVITY MODES
rez.time = obj.time;
rez.psth = obj.psth;
rez.condition = params.condition;
rez.alignEvent = params.alignEvent;

%% cd early

cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
epoch = 'early';
rez.cdEarly_mode = calcCD(obj,params,cond,epoch,rez.alignEvent);
clear cond


%% cd late mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
epoch = 'latedelay';
rez.cdLate_mode = calcCD(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% cd go mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
epoch = 'go';
rez.cdGo_mode = calcCD(obj,params,cond,epoch,rez.alignEvent);
clear cond


%% PLOTS

clrs = getColors();


% plot correct trials alone
plt.title = 'Correct Trials';
plt.legend = {'Right Hit','Left Hit'};
plt.conditions = [1,2];
plt.lw = [3 3];
plt.smooth = 31;
plt.colors = {clrs.rhit,clrs.lhit};
plt.save = 0;
plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 

% plot correct trials and error trials
plt.title = 'Correct and Error Trials';
plt.legend = {'Right Hit','Left Hit','Right Error', 'Left Error'};
plt.conditions = [1,2,5,6];
plt.lw = [2.5 2.5 1.5 1.5];
plt.smooth = 31;

plt.colors = {clrs.rhit,clrs.lhit,clrs.rmiss,clrs.lmiss};

plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 







