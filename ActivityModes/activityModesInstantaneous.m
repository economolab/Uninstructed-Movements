clear,clc,close all

if ispc
    pth = 'C:\Code\activity_modes';
elseif ismac
    pth = '/Users/Munib/Documents/Economo-Lab/code/ActivityModes';
end

addpath(genpath(pwd))

% had an idea and trying it out here
% find instantaneous mode at each time point rather than from the mean
% during an epoch
% then, do pca on that matrix of CVs. The top PCs are then the CVs. 

% didn't work very well :(

%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'goCue' or 'moveOnset'

params.lowFR               = 1; % remove clusters firing less than this val

% set conditions to use for projections
params.condition(1) = {'R&hit&~stim.enable&autowater.nums==2&~early'}; % right hits, no stim, aw off
params.condition(2) = {'L&hit&~stim.enable&autowater.nums==2&~early'}; % left hits, no stim, aw off
params.condition(3) = {'R&miss&~stim.enable&autowater.nums==2&~early'};   % error right, no stim, aw off
params.condition(4) = {'L&miss&~stim.enable&autowater.nums==2&~early'};   % error left, no stim, aw off
params.condition(5) = {'R&hit&~stim.enable&autowater.nums==1&~early'}; % right hits, no stim, aw on
params.condition(6) = {'L&hit&~stim.enable&autowater.nums==1&~early'}; % left hits, no stim, aw on
params.condition(7) = {'~hit&~miss&~stim.enable&autowater.nums==2&~early'}; % ignore

% set conditions used for finding the modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};

%% SET METADATA
% experiment meta data
if ispc
    meta.datapth = fullfile('C:\Code','data');
else
    meta.datapth = fullfile('/Users/Munib/Documents/Economo-Lab/code','data');
end
meta.anm = 'JEB7';
meta.date = '2021-04-29';
meta.datafn = findDataFn(meta);

meta.probe = 1;

% analysis meta data
meta.tmin = -3; % (s) relative to params.evName
meta.tmax = 3;  % (s) relative to params.evName
meta.dt = 0.005;

meta.smooth = 15; % smooth psth

% clusters (these qualities are included)
meta.quality = {'Fair','Good','Great','Excellent','single','multi'}; 

%% LOAD DATA
dat = load(fullfile(meta.datapth, meta.datafn));
obj = dat.obj;

obj.condition = params.condition;

%% get trials and clusters to use
meta.trialid = findTrials(obj, obj.condition);

cluQuality = {obj.clu{meta.probe}(:).quality}';
meta.cluid = findClusters(cluQuality, meta.quality);

%% align data
obj = alignSpikes(obj,meta,params);

%% get trial avg psth and single trial data
obj = getSeq(obj,meta);

%% remove low fr clusters
[obj, meta] = removeLowFRClusters(obj,meta,params);

%% ACTIVITY MODES
rez.time = obj.time;
rez.psth = obj.psth;
rez.condition = obj.condition;
rez.alignEvent = params.alignEvent;

%% choice mode (average over epoch)
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'delay';
rez.choice_mode = choiceMode(obj,meta,cond,epoch,rez.alignEvent);
clear cond

%% choice mode (instantaneous + pca)
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'delay';
rez.choice_mode_inst = choiceModeInst(obj,meta,cond,epoch,rez.alignEvent);
clear cond

% do some pca

%% orthogonalize

[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
modes = zeros(numel(meta.cluid),numel(fns));
for i = 1:numel(fns)
    modes(:,i) = rez.(fns{i});
end

orthModes = gschmidt(modes);

for i = 1:numel(fns)
    rez.(fns{i}) = orthModes(:,i);
end

%% PLOTS

% MODES VIZ

% plot correct trials alone
plt.title = 'Correct Trials';
plt.legend = {'Right Hit','Left Hit'};
plt.conditions = [1,2];
plt.lw = [2 2];
plt.smooth = 31;
plt.colors = {[0 0 1],[1 0 0]};
plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 

% plot correct trials and error trials
plt.title = 'Correct and Error Trials';
plt.legend = {'Right Hit','Left Hit','Right Error', 'Left Error'};
plt.conditions = [1,2,3,4];
plt.lw = [2.5 2.5 1.5 1.5];
plt.smooth = 31;
plt.colors = {[0 0 1],[1 0 0], ...
                 [0.5 0.5 1],[1 0.5 0.5]};
plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 

% plot correct trials and AW trials
plt.title = '2AFC and Autowater (Correct) Trials';
plt.legend = {'Right 2AFC','Left 2AFC','Right AW', 'Left AW'};
plt.conditions = [1,2,5,6];
plt.lw = [2.5 2.5 1.5 1.5];
plt.smooth = 31;
plt.colors = {[0 0 1],[1 0 0], ...
                 [0.2 0.8 0.9],[0.9 0.5 0.2]};
plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 

% plot correct and ignore trials
plt.title = 'Correct and Ignore Trials';
plt.legend = {'Right Hit','Left Hit','Ignore'};
plt.conditions = [1,2,7];
plt.lw = [2 2 2];
plt.smooth = 31;
plt.colors = {[0 0 1],[1 0 0],[0 0 0]};
plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 


% ORTHOGONALITY VIZ
% dotProductModes(rez,modes,'NOT ORTHOGONALIZED')
% dotProductModes(rez,orthModes,'ORTHOGONALIZED')

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fn = findDataFn(meta)
contents = dir(meta.datapth);
contents = {contents.name}';

strToFind = {'data_structure' , meta.anm, meta.date};

[fn,~] = patternMatchCellArray(contents, strToFind, 'all');
fn = fn{1};

end % loadRawDataObj












