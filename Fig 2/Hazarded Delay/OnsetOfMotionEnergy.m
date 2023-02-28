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
params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this mu

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
ctrlparams.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};        % All 2AFC hits, no stim

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
%% Load motion energy for control sessions
% % ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(ctrlmeta)
    disp(['Loading ME for control session ' num2str(sessix)])
    ctrlme(sessix) = loadMotionEnergy(ctrlobj(sessix), ctrlmeta(sessix), ctrlparams(sessix), datapth);
end
%% Load motion energy for hazarded delay sessions
camview = 'sidecam';
for sessix = 1:numel(meta)
    disp(['Loading ME for haz session ' num2str(sessix)])
    if obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{1},2) || obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{2},2)
        disp('Number of Bpod trials does not match number of video files for this session')
    else
        me(sessix) = loadMotionEnergy_Behav(obj(sessix), meta(sessix), params(sessix), datapth,camview);
    end
end
%% Load kinematic data
nSessions = numel(ctrlmeta);
for sessix = 1:numel(ctrlmeta)
    message = strcat('----Getting kinematic data for ctrl session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    ctrlkin(sessix) = getKinematics(ctrlobj(sessix), ctrlme(sessix), ctrlparams(sessix));
end
%% Load kin data for Haz animals
nSessions = numel(meta);
sess2incl = true(1,length(meta));
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for haz session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    if obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{1},2) || obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{2},2)
        disp('Number of Bpod trials does not match number of video files for this session')
        sess2incl(sessix) = 0;
    else
        kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
    end
end

% Remove sessions that aren't good
kin(~sess2incl) = [];
obj(~sess2incl) = [];
meta(~sess2incl) = [];
params(~sess2incl) = [];
%% Get rid of static delay sessions which don't have a delay period of 0.9
sess2incl = true(1,length(ctrlmeta));
for sessix = 1:numel(ctrlmeta)
    del = mode(ctrlobj(sessix).bp.ev.goCue)-mode(ctrlobj(sessix).bp.ev.delay);
    if del<0.85
        disp('Delay for this session is not 0.9')
        sess2incl(sessix) = 0;
    end
end
% Remove sessions that aren't good
ctrlkin(~sess2incl) = [];
ctrlobj(~sess2incl) = [];
ctrlmeta(~sess2incl) = [];
ctrlparams(~sess2incl) = [];
%% Get avg jaw velocity for all ctrl sessions
cond2use = 2;
feature = 'jaw_yvel_view1';


for sessix = 1:length(ctrlmeta)
    featix = find(strcmp(ctrlkin(sessix).featLeg,feature));
    condtrix = ctrlparams(sessix).trialid{cond2use};
    tempjaw = squeeze(ctrlkin(sessix).dat_std(:,condtrix,featix));
    tempjaw = abs(tempjaw);

%     imagesc(ctrlobj(1).time, 1:size(tempjaw,2),tempjaw')
%     pause

    ctrlobj(sessix).jawvel = mean(tempjaw,2);
end
%% Get avg jaw velocity for all haz delay sessions
feature = 'jaw_yvel_view1';
del2use = 1.2;

for sessix = 1:length(meta)
    featix = find(strcmp(kin(sessix).featLeg,feature));
    cond = obj(sessix).bp.hit&~obj(sessix).bp.stim.enable&~obj(sessix).bp.autowater&~obj(sessix).bp.early;
    delLength = obj(sessix).bp.ev.goCue-obj(sessix).bp.ev.delay;
    deltrix = find(delLength==del2use);
    condtrix = find(cond);
    touse = ismember(deltrix,condtrix);                     % Find which trials that are of the desired delay length are also in the desired condition
    ix = deltrix(touse);                                    % Get the trial numbers of trials which satisfy the condition above
    tempjaw = squeeze(kin(sessix).dat_std(:,ix,featix));
    tempjaw = abs(tempjaw);

    % Plot as a sanity check
%     imagesc(obj(1).time+0.5, 1:size(tempjaw,2),tempjaw')
%     pause
    obj(sessix).jawvel = mean(tempjaw,2);
end
%% Concatenate avg jaw velocities across all sessions
smooth = 21;
temp1 = [];
blah1 = [];

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
    temp1 = [temp1,mySmooth(ctrlobj(sessix).jawvel,smooth)];
    
end
jv.ctrl = temp1;

%%
delay = 0;
go = del2use;
startix = find(obj(1).time>delay,1,'first');
stopix = find(obj(1).time<go,1,'last');

temp2 = []; blah2 = [];
for sessix = 1:length(meta)
%     tempjaw = obj(sessix).jawvel;
%     deljaw = tempjaw(startix:stopix,:);
%     maxjaw = max(deljaw);
%     tempjaw = tempjaw./maxjaw;
%     temp2 = [temp2,mySmooth(tempjaw,smooth)];
    temp2 = [temp2,mySmooth(obj(sessix).jawvel,smooth)];
    
end
jv.haz = temp2;

clear temp1; clear temp2
%% Sanity check
% subplot(1,2,1)
% imagesc(ctrlobj(1).time,1:16,jv.ctrl')
% xlim([-1.3 2])
% colorbar
% 
% subplot(1,2,2)
% imagesc(obj(1).time+0.5,1:12,jv.haz')
% xlim([-1.3 2])
% colorbar
 %% Sanity check
% for i = 1:12
%     samp = mode(obj(i).bp.ev.sample)-mode(obj(i).bp.ev.delay);
%     trialstart = mode(obj(i).bp.ev.bitStart)-mode(obj(i).bp.ev.delay);
%     disp(num2str(trialstart))
% plot(obj(1).time+0.5,jv.haz(:,i)); hold on;
% xline(0,'k--')
% %xline(0.9,'k--')
% xline(1.2,'k--')
% xline(-1.3,'k--')
% hold off;
% pause
% end
%% Plot
colors = getColors_Updated();
ctrlcol = colors.afc;
hazcol = [0.5 0.5 0.5];
go = mode(ctrlobj(1).bp.ev.goCue)-mode(ctrlobj(1).bp.ev.(ctrlparams(1).alignEvent));
del2use = 1.2;
smooth = 70;
alph = 0.2;
figure();
ax = gca;
toplot = mySmooth(mean(jv.haz,2),smooth);
nSess = size(jv.haz,2);
err = 1.96*(mySmooth(std(jv.haz,0,2),smooth)./sqrt(nSess));
shadedErrorBar(obj(1).time+0.5,toplot,err,{'Color',hazcol,'LineWidth',2},alph,ax);
hold on;

toplot = mySmooth(mean(jv.ctrl,2),smooth);
nSess = size(jv.ctrl,2);
err = 1.96*(mySmooth(std(jv.ctrl,0,2),smooth)./sqrt(nSess));
shadedErrorBar(obj(1).time,toplot,err,{'Color',ctrlcol,'LineWidth',2},alph,ax);
xline(0,'LineStyle','--','Color','black','LineWidth',1.5)
xline(go,'LineStyle','--','Color',ctrlcol,'LineWidth',1.5)
xline(del2use,'LineStyle','--','Color',hazcol,'LineWidth',1.5)
legend(' ',' ',' ','Haz',' ',' ',' ','Static','Delay','Go cue, static','Go cue, haz','Location','best')
xlim([-1.3 1.5])
xlabel('Time from delay onset (s)')
ylabel([feature])
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