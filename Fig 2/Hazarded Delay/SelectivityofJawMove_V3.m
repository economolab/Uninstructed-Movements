% Compare jaw movement and velocity during the static delay and hazarded
% delay tasks
clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 3';
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Hazarded Delay')));
%% LOAD Hazarded Delay data
% Everything aligned to the delay period; tmin = -2.5, tmax = 3; dt = 1/200
% HAVE TO ACCOUNT FOR 0.5 SECOND DELAY WITHOUT SPIKEGLX
load('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 3\Hazarded Delay\BehavData\JEB11_JEB12_objs.mat');
load('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 3\Hazarded Delay\BehavData\JEB11_JEB12_meta.mat');
load('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 3\Hazarded Delay\BehavData\JEB11_JEB12_params.mat');

% Set params for hazarded delay data
params.tmin = meta(1).tmin; params.tmax = meta(1).tmax;
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw'}};
params.smooth = 15;
params.advance_movement = 0;

% Re-structure Hazarded Delay data to be the same format
[params,obj] = cleanUpData(params, meta,objs);
%% PARAMETERS
ctrlparams.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
ctrlparams.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
ctrlparams.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

ctrlparams.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
ctrlparams.condition(1)     = {'(hit|miss|no)'};                             % all trials
ctrlparams.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};        % R 2AFC hits, no stim
ctrlparams.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};        % L 2AFC hits, no stim

ctrlparams.tmin = -2.5;
ctrlparams.tmax = 3;
ctrlparams.dt = 1/200;

% smooth with causal gaussian kernel
ctrlparams.smooth = 15;

% cluster qualities to use
ctrlparams.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

ctrlparams.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

ctrlparams.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

ctrlparams.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

ctrlparams.advance_movement = 0;
%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

ctrlmeta = [];

% --- ALM ---
ctrlmeta = loadJEB6_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB7_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadEKH1_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadEKH3_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJGR2_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJGR3_ALMVideo(ctrlmeta,datapth);
%ctrlmeta = loadJEB13_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB14_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB15_ALMVideo(ctrlmeta,datapth);

ctrlparams.probe = {ctrlmeta.probe}; % put probe numbers into ctrlparams, one entry for element in ctrlmeta, just so i don't have to change code i've already written

%% LOAD DATA
% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% ctrlparams (struct array) - one entry per session
% ----------------------------------------------
[ctrlobj,ctrlparams] = loadSessionData_VidOnly(ctrlmeta,ctrlparams);
for sessix = 1:length(ctrlobj)
    ctrlobj(sessix).time = obj(1).time;
end
%%
% % ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(ctrlmeta)
    disp(['Loading ME for session ' num2str(sessix)])
    me(sessix) = loadMotionEnergy(ctrlobj(sessix), ctrlmeta(sessix), ctrlparams(sessix), datapth);
end
%% Load kinematic data
nSessions = numel(ctrlmeta);
for sessix = 1:numel(ctrlmeta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    ctrlkin(sessix) = getKinematics(ctrlobj(sessix), me(sessix), ctrlparams(sessix));
end
%% Load kin data for Haz animals
nSessions = numel(meta);
sess2incl = true(1,length(meta));
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    if obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{1},2) || obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{2},2)
        disp('Number of Bpod trials does not match number of video files for this session')
        sess2incl(sessix) = 0;
    else
        kin(sessix) = getKinematics_NoME(obj(sessix), params(sessix));
    end
end

% Remove sessions that aren't good
kin(~sess2incl) = [];
obj(~sess2incl) = [];
meta(~sess2incl) = [];
params(~sess2incl) = [];
%% Get avg jaw velocity for all ctrl sessions
cond2use = [2,3];
feature = 'jaw_yvel_view1';

for sessix = 1:length(ctrlmeta)
    featix = find(strcmp(ctrlkin(sessix).featLeg,feature));
    sessavg = NaN(length(ctrlobj(1).time),length(cond2use));
    for c = 1:length(cond2use)
        cond = cond2use(c);
        condtrix = ctrlparams(sessix).trialid{cond};
        tempjaw = squeeze(ctrlkin(sessix).dat(:,condtrix,featix));
        tempjaw = abs(tempjaw);
        sessavg(:,c) = mean(tempjaw,2);
    end
    ctrlobj(sessix).jawvel = sessavg;
end
%% Get avg jaw velocity for all haz delay sessions
feature = 'jaw_yvel_view1';
del2use = 1.2;

for sessix = 1:length(meta)
    featix = find(strcmp(kin(sessix).featLeg,feature));
    sessavg = NaN(length(ctrlobj(1).time),2);
    for c = 1:2
        if c==1
            cond = obj(sessix).bp.R&obj(sessix).bp.hit&~obj(sessix).bp.stim.enable&~obj(sessix).bp.autowater&~obj(sessix).bp.early;
        elseif c==2
            cond = obj(sessix).bp.L&obj(sessix).bp.hit&~obj(sessix).bp.stim.enable&~obj(sessix).bp.autowater&~obj(sessix).bp.early;
        end
        delLength = obj(sessix).bp.ev.goCue-obj(sessix).bp.ev.delay;
        deltrix = find(delLength==del2use);
        condtrix = find(cond);
        ix = find(ismember(deltrix,condtrix));
        tempjaw = squeeze(kin(sessix).dat(:,ix,featix));
        tempjaw = abs(tempjaw);
        sessavg(:,c) = mean(tempjaw,2);
    end
    obj(sessix).jawvel = sessavg;
end
%% Concatenate avg jaw velocities across all sessions
smooth = 41;
temp1 = [];

delay = mode(ctrlobj(1).bp.ev.delay)-mode(ctrlobj(1).bp.ev.(ctrlparams(1).alignEvent));
go = mode(ctrlobj(1).bp.ev.goCue)-mode(ctrlobj(1).bp.ev.(ctrlparams(1).alignEvent));
startix = find(ctrlobj(1).time>delay,1,'first');
stopix = find(ctrlobj(1).time<go,1,'last');
for sessix = 1:length(ctrlmeta)
    tempjaw = ctrlobj(sessix).jawvel;
     % Normalize to max jaw vel during delay period
%     deljaw = tempjaw(startix:stopix,:);
%     maxjaw = max(deljaw);
%     tempjaw = tempjaw./maxjaw;
%     temp1 = [temp1,mySmooth(tempjaw,smooth)];
    tempjaw = mySmooth(ctrlobj(sessix).jawvel(:,1),smooth)-mySmooth(ctrlobj(sessix).jawvel(:,2),smooth);
     if mean(tempjaw(startix:stopix))<0
        tempjaw = -1*tempjaw;
    end
    temp1 = [temp1,tempjaw];
end
jv.ctrl = temp1;
%%
delay = 0;
go = del2use;
startix = find(obj(1).time>delay,1,'first');
stopix = find(obj(1).time<go,1,'last');

temp2 = [];
for sessix = 1:length(meta)
%     tempjaw = obj(sessix).jawvel;
%     deljaw = tempjaw(startix:stopix,:);
%     maxjaw = max(deljaw);
%     tempjaw = tempjaw./maxjaw;
%     temp2 = [temp2,mySmooth(tempjaw,smooth)];
    %temp2 = [temp2,mySmooth(obj(sessix).jawvel,smooth)];
    tempjaw = mySmooth(obj(sessix).jawvel(:,1),smooth)-mySmooth(obj(sessix).jawvel(:,2),smooth);
    %plot(mySmooth(obj(sessix).jawvel(:,1),smooth),'blue'); hold on; plot(mySmooth(obj(sessix).jawvel(:,2),smooth)); hold off
    if mean(tempjaw(startix:stopix))<0
        tempjaw = -1*tempjaw;
    end
    temp2 = [temp2,tempjaw];
end
jv.haz = temp2;

clear temp1; clear temp2
%%
plot(obj(1).time+0.5,mean(jv.haz,2)); 
hold on; plot(obj(1).time,mean(jv.ctrl,2));
legend('Haz','Static')
xline(0.9)
xline(1.2)
xlim([-2.3 1.2])
%% FUNCTIONS
function [params,obj] = cleanUpData(params, meta,objs)
% Reorganize params into a struct to match structure of updated data
temp = params;
clear params
for sessix = 1:numel(meta)
    temp2 = temp;
    temp2.trialid = meta(sessix).trialid;
    params(sessix) = temp2;
end
clear temp temp2

% Reorganize objs into a struct to match structure of updated data
temp = objs;
clear objs
for sessix = 1:numel(meta)
    temp2 = temp{sessix};
    temp2.time = params(sessix).taxis; %+0.5
    obj(sessix) = temp2;
end
end