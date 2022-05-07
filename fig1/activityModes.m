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

% 0. context mode : defined during presample period 
%       (hit2afc - hitaw) / sqrt(sum(sd for each tt ^2));
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
% 7. remainder modes

%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 0.5; % remove clusters with firing rates across all trials less than this val

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

%% context mode
% hit2afc - hitaw during presample period
% cond{1} = params.condition{end-1}; % hit 2afc
% cond{2} = params.condition{end}; % hit aw
% rez.context_mode = (obj.presampleFR(:,end-1) - obj.presampleFR(:,end)) ./ sqrt(sum(obj.presampleSigma(:,(end-1):end).^2,2));
% clear cond

%% stimulus mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'sample'; 
rez.stimulus_mode = stimulusMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% choice mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'delay';
rez.choice_mode = choiceMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% action mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
epoch = 'action';
rez.action_mode = actionMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% outcome mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'outcome';
rez.outcome_mode = outcomeMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% ramping mode
cond{1} = params.modecondition{5};
epoch = {'presample','delay'};
rez.ramping_mode = rampingMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% go mode
cond{1} = params.modecondition{5};
epoch = {'postgo','prego'};
rez.go_mode = goMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% response mode 
% cond{1} = params.modecondition{1};
% cond{2} = params.modecondition{2};
% psthcond = [1,2];
% epoch = 'presample'; % used to estimate baseline firing rate
% rez.response_mode = responseMode(obj,params,cond,epoch,rez.alignEvent,psthcond);
% clear cond


%% orthogonalize

[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
modes = zeros(numel(params.cluid),numel(fns));
for i = 1:numel(fns)
    modes(:,i) = rez.(fns{i});
end

orthModes = gschmidt(modes);

for i = 1:numel(fns)
    rez.(fns{i}) = orthModes(:,i);
end

%% remainder mode
modesToKeep = eye(size(obj.psth,2)) - (orthModes*orthModes');

residualpsth = nan(size(obj.psth));
for i = 1:size(obj.psth,3)
    residualpsth(:,:,i) = obj.psth(:,:,i) * modesToKeep;
end

X = [residualpsth(:,:,1) ; residualpsth(:,:,2)]; % left and right 2afc

% SVD
V = myPCA(X - mean(X));
rez.remainder_mode = V(:,1); % S returned in decreasing order
% rez.remainder2_mode = V(:,2); % S returned in decreasing order
% rez.remainder3_mode = V(:,3); % S returned in decreasing order


%% PLOTS

% MODES VIZ

% plot correct trials alone
plt.title = 'Correct Trials';
plt.legend = {'Right Hit','Left Hit'};
plt.conditions = [1,2];
plt.lw = [3 3];
plt.smooth = 31;
plt.colors = {[0 0.4470 0.7410],[0.6350 0.0780 0.1840]};
plt.save = 0;
plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 

% % plot correct trials and error trials
% plt.title = 'Correct and Error Trials';
% plt.legend = {'Right Hit','Left Hit','Right Error', 'Left Error'};
% plt.conditions = [1,2,5,6];
% plt.lw = [2.5 2.5 1.5 1.5];
% plt.smooth = 31;
% plt.colors = {[0 0 1],[1 0 0], ...
%                  [0.5 0.5 1],[1 0.5 0.5]};
% plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 

% % plot correct trials and AW trials
% % plt.title used if plt.save=1
% if params.timeWarp && strcmpi(params.alignEvent,'gocue')
%     plt.title = '2AFC-AW (Correct) Trials [Time Warped | Go Cue]';
% elseif params.timeWarp && strcmpi(params.alignEvent,'firstlick')
%     plt.title = '2AFC-AW (Correct) Trials [Time Warped | First Lick]';
% elseif strcmpi(params.alignEvent,'gocue')
%     plt.title = '2AFC-AW (Correct) Trials [Go Cue]';
% elseif strcmpi(params.alignEvent,'firstlick')
%     plt.title = '2AFC-AW (Correct) Trials [First Lick]';
% end
% plt.title = '';
% plt.legend = {'Right 2AFC','Left 2AFC','Right AW', 'Left AW'};
% plt.conditions = [1,2,3,4];
% plt.lw = [2.7 2.7 2.7 2.7];
% plt.smooth = 31;
% plt.colors = {[0 0 1],[1 0 0], ...
%                  [190, 3, 252]./255,[252, 190, 3]./255};
% plt.save = 0;
% plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 

% % plot correct and ignore trials
% plt.title = 'Correct and Ignore Trials';
% plt.legend = {'Right Hit','Left Hit','Ignore'};
% plt.conditions = [1,2,7];
% plt.lw = [2 2 2];
% plt.smooth = 31;
% plt.colors = {[0 0 1],[1 0 0],[0 0 0]};
% plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 


% ORTHOGONALITY VIZ
dotProductModes(rez,modes,'NOT ORTHOGONALIZED')
dotProductModes(rez,orthModes,'ORTHOGONALIZED')
