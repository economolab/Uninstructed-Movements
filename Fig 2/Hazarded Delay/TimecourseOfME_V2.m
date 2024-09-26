% Compare jaw movement and velocity during the static delay and hazarded
% delay tasks

clear,clc,close all

whichcomp = 'LabPC';                                                % LabPC or Laptop

% Base path for code depending on laptop or lab PC
if strcmp(whichcomp,'LabPC')
    basepth = 'C:\Code';
elseif strcmp(whichcomp,'Laptop')
    basepth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code';
end


% add paths for data loading scripts, all fig funcs, and utils
utilspth = [basepth '\Munib Uninstruct Move\uninstructedMovements_v2'];
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
figpth = [basepth  '\Uninstructed-Movements\Fig 2'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Hazarded Delay')));
%% LOAD Hazarded Delay data
if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data\Hazarded Delay\BehavData';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 3\Hazarded Delay\BehavData';
end
    

% Everything aligned to the delay period; tmin = -2.5, tmax = 3; dt = 1/200
% HAVE TO ACCOUNT FOR 0.5 SECOND DELAY WITHOUT SPIKEGLX
load([datapth '\JEB11_JEB12_objs.mat']);
load([datapth '\JEB11_JEB12_meta.mat']);
load([datapth '\JEB11_JEB12_params.mat']);
%%
% Set params for hazarded delay data
params.tmin = meta(1).tmin; params.tmax = meta(1).tmax;
params.traj_features = {{'jaw'},...
    {'jaw'}};
params.smooth = 15;
params.advance_movement = 0;
params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this mu

% Re-structure Hazarded Delay data to be the same format
[params,obj] = cleanUpData(params, meta,objs);
%% Load motion energy for hazarded delay sessions
camview = 'sidecam';
whichcomp = 'LabPC';
if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

for sessix = 1:numel(meta)
    disp(['Loading ME for haz session ' num2str(sessix)])
    if obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{1},2) || obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{2},2)
        disp('Number of Bpod trials does not match number of video files for this session')
    else
        me(sessix) = loadMotionEnergy_Behav(obj(sessix), meta(sessix), params(sessix), datapth,camview);
    end
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
%% Get avg motion energy for all haz delay sessions
feature = 'motion_energy';
del2use = 1.2;
sm = 100;

startix = find(obj(1).time>0,1,'first');
stopix = find(obj(1).time<del2use,1,'last');

for sessix = 1:length(meta)
    featix = find(strcmp(kin(sessix).featLeg,feature));
    tempjaw = [];
    for c = 1:2
        if c==1
            cond = obj(sessix).bp.R&obj(sessix).bp.hit&~obj(sessix).bp.stim.enable&~obj(sessix).bp.autowater&~obj(sessix).bp.early;
        else
            cond = obj(sessix).bp.L&obj(sessix).bp.hit&~obj(sessix).bp.stim.enable&~obj(sessix).bp.autowater&~obj(sessix).bp.early;

        end
        delLength = obj(sessix).bp.ev.goCue-obj(sessix).bp.ev.delay;
        deltrix = find(delLength==del2use);
        condtrix = find(cond);
        touse = ismember(deltrix,condtrix);                     % Find which trials that are of the desired delay length are also in the desired condition
        ix = deltrix(touse);                                    % Get the trial numbers of trials which satisfy the condition above
        temp = squeeze(kin(sessix).dat(:,ix,featix));
        temp = abs(temp);
        condtemp = mySmooth(mean(temp,2,'omitnan'),41);
        tempjaw = [tempjaw, condtemp];
    end

    sel = tempjaw(:,1)-tempjaw(:,2);
    delSelectivity = mean(sel(startix:stopix),'omitnan');
    if delSelectivity<0
        sel = -1*sel;
    end
    obj(sessix).jawSelectivity = mySmooth(sel,sm);
    obj(sessix).jawvel = mean(tempjaw,2);
end
%% Concatenate avg motion energy across all rand sessions
delay = 0;
go = del2use;
startix = find(obj(1).time>delay,1,'first');
stopix = find(obj(1).time<go,1,'last');

temp2 = []; blah2 = [];
for sessix = 1:length(meta)
%     temp2 = [temp2,mySmooth(obj(sessix).jawSelectivity,smooth)];
    temp2 = [temp2,obj(sessix).jawvel];
    
end
jv.haz = temp2;

clear temp1; clear temp2
clearvars -except jv
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------- JEB23 and JEB24 data (different format than JEB11 and JEB12) --------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS 
paramsHaz.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
paramsHaz.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
paramsHaz.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

paramsHaz.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
paramsHaz.condition(1)     = {'(hit|miss|no)'};                             % all trials
paramsHaz.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};      % All 2AFC hits, no stim
paramsHaz.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};
paramsHaz.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};

paramsHaz.tmin = -2.5;
paramsHaz.tmax = 3;
paramsHaz.dt = 1/200;

% smooth with causal gaussian kernel
paramsHaz.smooth = 15;

% cluster qualities to use
paramsHaz.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

paramsHaz.traj_features = {{'jaw'},...
    {'jaw'}};

paramsHaz.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

paramsHaz.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

paramsHaz.advance_movement = 0;
%% SPECIFY DATA TO LOAD
whichcomp = 'LabPC';

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

% Base path for code depending on laptop or lab PC
if strcmp(whichcomp,'LabPC')
    basepth = 'C:\Code';
elseif strcmp(whichcomp,'Laptop')
    basepth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code';
end


metaHaz = [];

metaHaz = loadJEB23_ALMVideo(metaHaz,datapth);
metaHaz = loadJEB24_ALMVideo(metaHaz,datapth);

paramsHaz.probe = {metaHaz.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

%% LOAD DATA
% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[objHaz,paramsHaz] = loadSessionData(metaHaz,paramsHaz);
%%
% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(metaHaz)
    disp(['Loading ME for session ' num2str(sessix)])
    meHaz(sessix) = loadMotionEnergy(objHaz(sessix), metaHaz(sessix), paramsHaz(sessix), datapth);
end
%% Load kinematic data
nSessions = numel(metaHaz);
for sessix = 1:numel(metaHaz)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kinHaz(sessix) = getKinematics(objHaz(sessix), meHaz(sessix), paramsHaz(sessix));
end
%%
del2use = 1.2000;
sm = 50;
cond2use = 4;

delay = -1.3;
start = -1.5;
startix = find(objHaz(1).time>start,1,'first');
stopix = find(objHaz(1).time<delay,1,'last');

kinfeat = 'motion_energy';
featix = find(strcmp(kinHaz(1).featLeg,kinfeat));

for sessix = 1:length(objHaz)
    condtrix = paramsHaz(sessix).trialid{cond2use};
    delLength = objHaz(sessix).bp.ev.goCue-objHaz(sessix).bp.ev.delay;
    deltrix = find(delLength<1.3&delLength>1.1);
    touse = ismember(deltrix,condtrix);                     % Find which trials that are of the desired delay length are also in the desired condition
    ix = deltrix(touse);                                    % Get the trial numbers of trials which satisfy the condition above
    delkin = kinHaz(sessix).dat(:,ix,featix);
    meankin = squeeze(mean(delkin,2,'omitnan'));
    baseME = mean(meankin(startix:stopix,:),1,'omitnan');
    meankin = meankin-baseME;
    randME(:,sessix) = mySmooth(meankin,sm);
end
%% Group session-averaged motion energy values from JEB11,12 and JEB23,24

%%%%%%% HAVE TO SHIFT JEB11 and 12 by 0.5 s before combining %%%%%%%%%%%%%
time1 = -1.5-0.5;
time2 = 1.2-0.5;
ix1 = find(objHaz(1).time>time1,1,'first');
ix2 = find(objHaz(1).time<time2,1,'last');

time3 = -1.5;
time4 = 1.2;
ix3 = find(objHaz(1).time>time3,1,'first');
ix4 = find(objHaz(1).time<time4,1,'last');

% plot(objHaz(1).time(ix3:ix4), mean(jv.haz(ix1:ix2,:),2,'omitnan')); hold on; plot(objHaz(1).time(ix3:ix4), mean(randME(ix3:ix4,:),2,'omitnan'));
% legend('JEBs11/12','JEBs23/24')

JEB1112_jv = jv.haz;
jv.haz = [jv.haz(ix1:ix2,:),randME(ix3:ix4,:)];

HazTimeIx1 = ix3;
HazTimeIx2 = ix4;

clearvars -except jv randME JEB1112_jv objHaz HazTimeIx1 HazTimeIx2
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------FIXED DELAY SESSIONS---------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
ctrlparams.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
ctrlparams.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
ctrlparams.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

ctrlparams.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
ctrlparams.condition(1)     = {'(hit|miss|no)'};                             % all trials
ctrlparams.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};        % All 2AFC hits, no stim
ctrlparams.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};

ctrlparams.tmin = -2.5;
ctrlparams.tmax = 3;
ctrlparams.dt = 1/200;

% smooth with causal gaussian kernel
ctrlparams.smooth = 15;

% cluster qualities to use
ctrlparams.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

ctrlparams.traj_features = {{'jaw'},...
    {'jaw'}};

ctrlparams.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

ctrlparams.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

ctrlparams.advance_movement = 0;
%% SPECIFY DATA TO LOAD
whichcomp = 'LabPC';

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end


ctrlmeta = [];

% --- ALM ---
ctrlmeta = loadJEB6_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB7_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadEKH1_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadEKH3_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJGR2_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJGR3_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB13_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB14_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB15_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB19_ALMVideo(ctrlmeta,datapth);

ctrlparams.probe = {ctrlmeta.probe}; % put probe numbers into ctrlparams, one entry for element in ctrlmeta, just so i don't have to change code i've already written
%% LOAD DATA
% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% ctrlparams (struct array) - one entry per session
% ----------------------------------------------
[ctrlobj,ctrlparams] = loadSessionData_VidOnly(ctrlmeta,ctrlparams);
for sessix = 1:length(ctrlobj)
    ctrlobj(sessix).time = objHaz(1).time;
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

%% Load kinematic data
nSessions = numel(ctrlmeta);
for sessix = 1:numel(ctrlmeta)
    message = strcat('----Getting kinematic data for ctrl session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    ctrlkin(sessix) = getKinematics(ctrlobj(sessix), ctrlme(sessix), ctrlparams(sessix));
end
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
%% Get avg motion energy for all ctrl sessions
cond2use = [2,3];
feature = 'motion_energy';
sm = 100;

samp = mode(ctrlobj(1).bp.ev.sample)-mode(ctrlobj(1).bp.ev.(ctrlparams(1).alignEvent));
delay = mode(ctrlobj(1).bp.ev.delay)-mode(ctrlobj(1).bp.ev.(ctrlparams(1).alignEvent));
go = mode(ctrlobj(1).bp.ev.goCue)-mode(ctrlobj(1).bp.ev.(ctrlparams(1).alignEvent));
startix = find(ctrlobj(1).time>delay,1,'first');
stopix = find(ctrlobj(1).time<go,1,'last');

for sessix = 1:length(ctrlmeta)
    featix = find(strcmp(ctrlkin(sessix).featLeg,feature));
    tempjaw = [];
    for c = 1:length(cond2use)
        cond = cond2use(c);
        condtrix = ctrlparams(sessix).trialid{cond};
        temp = squeeze(ctrlkin(sessix).dat(:,condtrix,featix));
        temp = abs(temp);
        condtemp = mySmooth(mean(temp,2,'omitnan'),sm);
        tempjaw = [tempjaw, condtemp];
    end

    ctrlobj(sessix).jawvel = mean(tempjaw,2);
    sel = tempjaw(:,1)-tempjaw(:,2);
    delSelectivity = mean(sel(startix:stopix),'omitnan');
    if delSelectivity<0
        sel = -1*sel;
    end
    ctrlobj(sessix).jawSelectivity = mySmooth(sel,sm);
end
%% Concatenate avg motion energy across all ctrl sessions
%%% Have to account for the fact that there is a 0.5 s shift in data from JEB19
time1 = -1.5+0.5;
time2 = 1.2+0.5;
ix1 = find(objHaz(1).time>time1,1,'first');
ix2 = find(objHaz(1).time<time2,1,'last');

time3 = -1.5;
time4 = 1.2;
ix3 = find(objHaz(1).time>time3,1,'first');
ix4 = find(objHaz(1).time<time4,1,'last');

temp1 = [];
for sessix = 1:length(ctrlmeta)
    anmName = ctrlmeta(sessix).anm;
    % Shifts JEB19 data back by 0.5 s
    if strcmp(anmName,'JEB19')
        tempjaw = ctrlobj(sessix).jawvel(ix1:ix2);
    else
        tempjaw = ctrlobj(sessix).jawvel(ix3:ix4);
    end
    temp1 = [temp1,tempjaw]; 
end
jv.ctrl = temp1;

%% Baseline subtract  Motion energy values to presample
updatedCtrlTime = ctrlobj(1).time(ix3:ix4);
start = -1.5;
stop = -1.3;
startix = find(updatedCtrlTime>start,1,'first');
stopix = find(updatedCtrlTime<stop,1,'last');

baseME = mean(jv.ctrl(startix:stopix,:),1,'omitnan');
baseME = mean(baseME,'omitnan');
jv.ctrl = jv.ctrl-baseME;
%% Normalize to the 99th percentile within delay to account for differences in ME magnitudes?
normval = 99;
p99haz = prctile(jv.haz,normval);
jv.haznorm = jv.haz./p99haz;

start = 0;
stop = 0.9;
startix = find(updatedCtrlTime>start,1,'first');
stopix = find(updatedCtrlTime<stop,1,'last');
p99fix = prctile(jv.ctrl(startix:stopix,:),normval);
jv.ctrlnorm = jv.ctrl./p99fix;
%% Compare the slopes of the two curves (not normalized)
smooth = 70;

nTimePts = size(jv.ctrl,1)-1;
nSess = size(jv.ctrl,2);
slope.fixed = NaN(nTimePts,nSess);
for sessix = 1:nSess
    tempslope = diff(jv.ctrl(:,sessix));
    slope.fixed(:,sessix) = mySmooth(tempslope,smooth);
end

nTimePts = size(jv.haz,1)-1;
nSess = size(jv.haz,2);
slope.rand = NaN(nTimePts,nSess);
for sessix = 1:nSess
    tempslope = diff(jv.haz(:,sessix));
    slope.rand(:,sessix) = mySmooth(tempslope,smooth);
end

%% Two-sided t-test to determine whether the slopes of the ME curves are different from one another at each time point
siglevel = 0.05;
nWindows = size(slope.rand,1);   % # of time windows being used for the t-test
zerodif = NaN(nWindows,1);            % [time x 1] store whether ttest2 accepts (0) or rejects (1) the null hyp at this time point
for tt = 1:nWindows
    fixedSel = slope.fixed(tt,:);
    randSel = slope.rand(tt,:);
    hyp = ttest2(fixedSel,randSel,"Alpha",siglevel);
    if hyp == 0
        hyp = NaN;
    end
    zerodif(tt)=hyp;
end
%% Two-sided t-test to determine whether the ME curves are different from one another at each time point
siglevel = 0.05;
nWindows = size(jv.ctrl,1);   % # of time windows being used for the t-test
zerodif = NaN(nWindows,1);            % [time x 1] store whether ttest2 accepts (0) or rejects (1) the null hyp at this time point
for tt = 1:nWindows
    fixedME = jv.ctrl(tt,:);
    randME = jv.haz(tt,:);
    hyp = ttest2(fixedME,randME,"Alpha",siglevel);
    if hyp == 0
        hyp = NaN;
    end
    zerodif(tt)=hyp;
end
%% Plot
clearvars -except ctrlkin ctrlmeta ctrlobj ctrlparams jv kin meta obj params feature ...
    jv randME JEB1112_jv objHaz HazTimeIx1 HazTimeIx2 updatedCtrlTime zerodif

colors = getColors();
ctrlcol = colors.afc;
hazcol = [0.5 0.5 0.5];
go = mode(ctrlobj(1).bp.ev.goCue)-mode(ctrlobj(1).bp.ev.(ctrlparams(1).alignEvent));
ctrlstop = find(updatedCtrlTime<go,1,'last');
del2use = 1.2;
smooth = 10;
alph = 0.2;

figure();
ax = gca;
toplot = mean(mySmooth(jv.haz,smooth),2);
nSess = size(jv.haz,2);
%err = 1.96*(mySmooth(std(jv.haz,0,2),smooth)./sqrt(nSess));
err = std(mySmooth(jv.haz,smooth),0,2)./sqrt(nSess);
shadedErrorBar(objHaz(1).time(HazTimeIx1:HazTimeIx2),toplot,err,{'Color',hazcol,'LineWidth',2},alph,ax);
hold on;

toplot = mean(mySmooth(jv.ctrl,smooth),2);
nSess = size(jv.ctrl,2);
% err = 1.96*(mySmooth(std(jv.ctrl,0,2),smooth)./sqrt(nSess));
err = std(mySmooth(jv.ctrl,smooth),0,2)./sqrt(nSess);
shadedErrorBar(updatedCtrlTime(1:ctrlstop),toplot(1:ctrlstop),err(1:ctrlstop),{'Color',ctrlcol,'LineWidth',2},alph,ax);

plot(updatedCtrlTime,45*zerodif(:,1),'Color',[0 0 0],'LineWidth',3)

xline(0,'LineStyle','--','Color','black','LineWidth',1.5)
xline(go,'LineStyle','--','Color',ctrlcol,'LineWidth',1.5)
xline(del2use,'LineStyle','--','Color',hazcol,'LineWidth',1.5)
legend(' ',' ',' ','Haz',' ',' ',' ','Static','Delay','Go cue, static','Go cue, haz','Location','best')
xlim([-1.3 del2use])
ylim([0 47])
xlabel('Time from delay onset (s)')
ylabel(['Average ' feature])
%%
colors = getColors();
ctrlcol = colors.afc;
hazcol = [0.5 0.5 0.5];
go = mode(ctrlobj(1).bp.ev.goCue)-mode(ctrlobj(1).bp.ev.(ctrlparams(1).alignEvent));
ctrlstop = find(updatedCtrlTime<go,1,'last');
del2use = 1.2;
smooth = 10;
alph = 0.2;

figure();
ax = gca;
toplot = mean(mySmooth(jv.haznorm,smooth),2);
nSess = size(jv.haz,2);
%err = 1.96*(mySmooth(std(jv.haz,0,2),smooth)./sqrt(nSess));
err = std(mySmooth(jv.haznorm,smooth),0,2)./sqrt(nSess);
shadedErrorBar(objHaz(1).time(HazTimeIx1:HazTimeIx2),toplot,err,{'Color',hazcol,'LineWidth',2},alph,ax);
hold on;

toplot = mean(mySmooth(jv.ctrlnorm,smooth),2);
nSess = size(jv.ctrlnorm,2);
% err = 1.96*(mySmooth(std(jv.ctrl,0,2),smooth)./sqrt(nSess));
err = std(mySmooth(jv.ctrlnorm,smooth),0,2)./sqrt(nSess);
shadedErrorBar(updatedCtrlTime(1:ctrlstop),toplot(1:ctrlstop),err(1:ctrlstop),{'Color',ctrlcol,'LineWidth',2},alph,ax);
xline(0,'LineStyle','--','Color','black','LineWidth',1.5)
xline(go,'LineStyle','--','Color',ctrlcol,'LineWidth',1.5)
xline(del2use,'LineStyle','--','Color',hazcol,'LineWidth',1.5)
legend(' ',' ',' ','Haz',' ',' ',' ','Static','Delay','Go cue, static','Go cue, haz','Location','best')
xlim([-1.3 del2use])
xlabel('Time from delay onset (s)')
ylabel(['Average ' feature])
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