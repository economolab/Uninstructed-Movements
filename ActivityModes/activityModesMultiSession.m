clear,clc,close all

if ispc
    pth = 'C:\Code\activity_modes';
elseif ismac
    pth = '/Users/Munib/Documents/Economo-Lab/code/activity_modes';
end

addpath(genpath(pwd))

% find all 7 activity modes as described in :
% Thalamus drives diverse responses in frontal cortex during decision-making
% Weiguo Yang, Sri Laasya Tipparaju, Guang Chen, Nuo Li

% 1. stimulus mode: defined during stimulus (sample) period
%       ((hitR - missL) + (missR - hitL)) / sqrt(sum(sd for each tt ^2));
% 2. choice mode: defined during delay period
%       ((hitR - missR) + (missL - hitL)) / sqrt(sum(sd for each tt ^2));
% 3. action mode: defined during mvmt init (0.1 to 0.3 s rel to go cue)
%       (hitR - hitL) / sqrt(sum(sd for each tt ^2));
% 4. outcome mode: defined during response epoch (0 to 1.3 s rel go cue)
%       ((hitR - missL) + (missR - hitL)) / sqrt(sum(sd for each tt ^2));
% 5. ramping mode: in correct trials
%       (hit_presample - hit_delay) / sqrt(sum(sd for each tt ^2));
% 6. go mode: 0.1 sec before or after go cue
%       (hit_after - hit_before) / sqrt(sum(sd for each tt ^2));
% 7. response mode: 
%    a. find eigenvectors of basline subtracted PSTHs using SVD
%       aa. matrix was of size (n x (2t)), where left and right trials concatenated
%       in time
%    b. response mode = eigenvector with largest eigenvalue


%% TODO
% subtract out modes found using 2afc from aw context. See what's left.
% preprocess data other than normalize???
% handle multiple probes

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
meta(1).datapth = fullfile('C:\Code','data');
meta(1).anm = 'JEB7';
meta(1).date = '2021-04-29';
meta(1).datafn = findDataFn(meta(1));
meta(1).probe = 1;
% analysis meta data
meta(1).tmin = -3; % (s) relative to params.alignEvent
meta(1).tmax = 3;  % (s) relative to params.alignEvent
meta(1).dt = 0.005;
meta(1).smooth = 15; % smooth psth
% clusters (these qualities are included)
meta(1).quality = {'Fair','Good','Great','Excellent','single','multi'}; 

meta(2) = meta(1); % use most of the same fields
meta(2).datapth = 'Y:\JEB\Experiments\JEB7\Analysis\2021-04-18';
meta(2).anm = 'JEB7';
meta(2).date = '2021-04-18';
meta(2).datafn = findDataFn(meta(2));


%% PREPROCESS DATA
objs = loadObjs(meta);

%% PREPROCESS DATA
for i = 1:numel(meta)
    obj = objs{i};
    obj.condition = params.condition;
    % get trials and clusters to use
    meta(i).trialid = findTrials(obj, obj.condition);
    cluQuality = {obj.clu{meta(i).probe}(:).quality}';
    meta(i).cluid = findClusters(cluQuality, meta(i).quality);
    % align data
    obj = alignSpikes(obj,meta(i),params);
    % get trial avg psth and single trial data
    obj = getSeq(obj,meta(i));
    % remove low fr clusters
    [obj, meta(i)] = removeLowFRClusters(obj,meta(i),params);
    objs{i} = obj;
end

%% ACTIVITY MODES
rez.time = objs{1}.time;
rez.psth = concatPSTH(objs);
rez.condition = objs{1}.condition;
rez.alignEvent = params.alignEvent;

%% stimulus mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'sample';
rez.stimulus_mode = stimulusModeMulti(objs,meta,cond,epoch,rez.alignEvent);
clear cond

%% choice mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'delay';
rez.choice_mode = choiceModeMulti(objs,meta,cond,epoch,rez.alignEvent);
clear cond

%% action mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
epoch = 'action';
rez.action_mode = actionModeMulti(objs,meta,cond,epoch,rez.alignEvent);
clear cond

%% outcome mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'outcome';
rez.outcome_mode = outcomeModeMulti(objs,meta,cond,epoch,rez.alignEvent);
clear cond

%% ramping mode
cond{1} = params.modecondition{5};
epoch = {'presample','delay'};
rez.ramping_mode = rampingModeMulti(objs,meta,cond,epoch,rez.alignEvent);
clear cond

%% go mode
cond{1} = params.modecondition{5};
epoch = {'postgo','prego'};
rez.go_mode = goModeMulti(objs,meta,cond,epoch,rez.alignEvent);
clear cond

%% response mode 
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
psthcond = [1,2];
epoch = 'presample'; % used to estimate baseline firing rate
rez.response_mode = responseModeMulti(objs,meta,cond,epoch,rez.alignEvent,rez.psth,psthcond);
clear cond

%% orthogonalize

[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
modes = zeros(size(rez.psth,2),numel(fns));
for i = 1:numel(fns)
    modes(:,i) = rez.(fns{i});
end

orthModes = gschmidt(modes);

for i = 1:numel(fns)
    rez.(fns{i}) = orthModes(:,i);
end

%% PLOTS

% LATENTS VIZ

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
function objfn = findDataFn(meta)
contents = dir(meta.datapth);
contents = {contents.name}';

strToFind = {'data_structure' , meta.anm, meta.date};

[fn,~] = patternMatchCellArray(contents, strToFind, 'all');
objfn = fn{1};

end % loadRawDataObj

function objs = loadObjs(meta)
objs = cell(1,numel(meta));
for i = 1:numel(meta)
    dat = load(fullfile(meta(i).datapth, meta(i).datafn));
    objs{i} = dat.obj;
end
end % loadObjs

function psth = concatPSTH(objs)
psth = objs{1}.psth;
for i = 2:numel(objs)
    psth = cat(2,psth,objs{i}.psth);
end
end % concatPSTH












