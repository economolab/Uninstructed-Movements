% Compare motion energy during the static delay and hazarded delay tasks
% (within animal comparison)

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
%% PARAMETERS
params.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};      % All 2AFC hits, no stim
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};

params.tmin = -2.5;
params.tmax = 3;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

params.traj_features = {{'jaw','tongue'},...
    {'jaw','tongue'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;
%% SPECIFY DATA TO LOAD

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end


meta = [];

meta = loadJEB23_Behav(meta,datapth);
meta = loadJEB24_Behav(meta,datapth);

%params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

%% LOAD DATA

[obj,params] = loadBehavSessionData(meta,params);
for sessix = 1:length(obj)
    disp(['Loading data object for session ' num2str(sessix)])
    obj(sessix).time = params(sessix).tmin:params(sessix).dt:params(sessix).tmax;
end
%%
% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    disp(['Loading ME for session ' num2str(sessix)])
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%%
% Ideally want to have 30 [numDesiredTrials] L and 30 R trials before infusion and 
% 30 L and 30 R trials after infusion

% Right now can get away with 15 of each since this is preliminary analysis

numDesiredTrials = 30;
for sessix = 1:length(meta)
    check(1) = length(params(sessix).trialid{2})>numDesiredTrials;      % Greater than [numDesired] L trials before infusion
    check(end+1) = length(params(sessix).trialid{3})>numDesiredTrials;  % Greater than [numDesired] R trials before infusion
    
    if sum(check) == length(check)
        disp(['Session ' meta(sessix).anm ' ' meta(sessix).date ' fits trial inclusion criteria'])
    else
        disp(['Session ' meta(sessix).anm ' ' meta(sessix).date ' DOES NOT FIT inclusion criteria'])
    end
end
clearvars -except meta obj params kin me
%% Get average ME for fixed delay sessions
sess2use = NaN(1,length(meta));
for ss = 1:length(meta)
    deltype = meta(ss).delType;
    if strcmp(deltype,'fixed')
        sess2use(ss) = 1;
    end
end
fixedsessix = find(sess2use==1);

avgME = NaN(length(obj(1).time),length(fixedsessix));
sm = 50;
cond2use = 4;
kinfeat = 'motion_energy';
featix = find(strcmp(kin(1).featLeg,kinfeat));

%%%%%%% Probably want this to be 99th percentile instead %%%%%%%%%%%%%%%%%%
% Baseline subtract so that ME is centered at zero (subtracting presample
% motion energy) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delay = -1.3;
start = -1.5;
% delay = 2.5;
% start = 1.5;
startix = find(obj(1).time>start,1,'first');
stopix = find(obj(1).time<delay,1,'last');

for ss = 1:length(fixedsessix)
    sessix = fixedsessix(ss);
    condtrix = params(sessix).trialid{cond2use};
    condkin = kin(sessix).dat(:,condtrix,featix);
    meankin = squeeze(mean(condkin,2,'omitnan'));
    baseME = mean(meankin(startix:stopix,:),1,'omitnan');
    meankin = meankin-baseME;
    avgME(:,ss) = mySmooth(meankin,sm);
end
%%
% Get session IDs that do not have fixed delay period 
sess2use = NaN(1,length(meta));
randDay = cell(1,length(meta));
for ss = 1:length(meta)
    deltype = meta(ss).delType;
    if ~strcmp(deltype,'fixed')
        sess2use(ss) = 1;
        randDay{ss} = deltype;
    end
   
end
randsessix = find(sess2use==1);

feature = 'motion_energy';

del2use = 1.2000;
sm = 50;
cond2use = 4;

delay = -1.3;
start = -1.5;
% delay = 1.3;
% start = 1.0;
startix = find(obj(1).time>start,1,'first');
stopix = find(obj(1).time<delay,1,'last');

for ss = 1:length(randsessix)
    sessix = randsessix(ss);
    condtrix = params(sessix).trialid{cond2use};
    delLength = obj(sessix).bp.ev.goCue-obj(sessix).bp.ev.delay;
    deltrix = find(delLength<1.3&delLength>1.1);
    touse = ismember(deltrix,condtrix);                     % Find which trials that are of the desired delay length are also in the desired condition
    ix = deltrix(touse);                                    % Get the trial numbers of trials which satisfy the condition above
    delkin = kin(sessix).dat(:,ix,featix);
    meankin = squeeze(mean(delkin,2,'omitnan'));
    baseME = mean(meankin(startix:stopix,:),1,'omitnan');
    meankin = meankin-baseME;
    randME(:,ss) = mySmooth(meankin,sm);

    randsessname{ss} = [meta(sessix).anm ';' meta(sessix).date];
end
%% 
plot(obj(1).time,mean(avgME,2,'omitnan'),'LineWidth',2,'Color',[0.5 0.5 0.5]); hold on;
plot(obj(1).time,mean(randME,2,'omitnan'),'LineWidth',2,'Color',[0.8 0.1 0.8]);
xlim([-1.4 2.5])
xline(0.9,'--','Color',[0.5 0.5 0.5],'LineWidth',1.5)
xline(1.2,'m--','LineWidth',1.5)
xline(0,'k-.','LineWidth',1.5)
%xline(-1.3,'g--','LineWidth',1.5)
xlabel('Time from delay onset (s)')
ylabel('Motion energy')
legend({'Fixed delay','Rand delay','Go cue - fixed','Go cue - randomized'},'Location','best')
%%
blah = {};
for rr = 1:length(randDay)
    if ~isempty(randDay{rr})
        blah = [blah, randDay{rr}];
    end
end
%%
cols = {[1 0 1],[1 0.25 1],[0.75 0.5 0.75]};
for ii = 1:(size(randME,2))
    cc = ii*0.1;
    col = cc*[1 0 1];
    plot(obj(1).time,randME(:,ii),'LineWidth',1.5,'Color',col); hold on
end
legend(randsessname)

plot(obj(1).time,mean(avgME,2,'omitnan'),'LineWidth',3,'Color',[0.6 0.6 0.6]); hold on;
xlim([-1.5 2.5])
xline(0,'k--','LineWidth',1.5)
%xline(-1.3,'g--','LineWidth',1.5)
xline(0.9,'r--','LineWidth',1.5)

xlabel('Time from delay onset (s)')
ylabel('Motion energy')
% legend(randDay)
%%
sessix = 1;
randix = find(randsessix==sessix);
sessname = randsessname{randix};
feat = 'motion_energy';


condtrix = params(sessix).trialid{cond2use};
delLength = obj(sessix).bp.ev.goCue-obj(sessix).bp.ev.delay;
deltrix = find(delLength<1.3&delLength>1.1);
touse = ismember(deltrix,condtrix);                     % Find which trials that are of the desired delay length are also in the desired condition
ix = deltrix(touse);                                    % Get the trial numbers of trials which satisfy the condition above
delkin = kin(sessix).dat(:,ix,featix);
nTrials = size(delkin,2);

figure()
imagesc(obj(1).time,1:nTrials,delkin')
xlim([-1.5 2.5])
title(sessname)
colorbar
%% Get avg motion energy for all ctrl sessions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cond2use = [2,3];
% feature = 'motion_energy';
% sm = 100;
% 
% samp = mode(obj(1).bp.ev.sample)-mode(obj(1).bp.ev.(params(1).alignEvent));
% delay = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent));
% go = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent));
% startix = find(obj(1).time>delay,1,'first');
% stopix = find(obj(1).time<go,1,'last');
% 
% for sessix = 1:length(meta)
%     featix = find(strcmp(ctrlkin(sessix).featLeg,feature));
%     tempjaw = [];
%     for c = 1:length(cond2use)
%         cond = cond2use(c);
%         condtrix = params(sessix).trialid{cond};
%         temp = squeeze(ctrlkin(sessix).dat(:,condtrix,featix));
%         temp = abs(temp);
%         condtemp = mySmooth(mean(temp,2,'omitnan'),sm);
%         tempjaw = [tempjaw, condtemp];
%     end
% 
%     obj(sessix).jawvel = mean(tempjaw,2);
%     sel = tempjaw(:,1)-tempjaw(:,2);
%     delSelectivity = mean(sel(startix:stopix),'omitnan');
%     if delSelectivity<0
%         sel = -1*sel;
%     end
%     obj(sessix).jawSelectivity = mySmooth(sel,sm);
% end
% %% Get avg jaw velocity for all haz delay sessions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% feature = 'motion_energy';
% del2use = 1.2;
% sm = 100;
% 
% startix = find(obj(1).time>0,1,'first');
% stopix = find(obj(1).time<del2use,1,'last');
% 
% for sessix = 1:length(meta)
%     featix = find(strcmp(kin(sessix).featLeg,feature));
%     tempjaw = [];
%     for c = 1:2
%         if c==1
%             cond = obj(sessix).bp.R&obj(sessix).bp.hit&~obj(sessix).bp.stim.enable&~obj(sessix).bp.autowater&~obj(sessix).bp.early;
%         else
%             cond = obj(sessix).bp.L&obj(sessix).bp.hit&~obj(sessix).bp.stim.enable&~obj(sessix).bp.autowater&~obj(sessix).bp.early;
% 
%         end
%         delLength = obj(sessix).bp.ev.goCue-obj(sessix).bp.ev.delay;
%         deltrix = find(delLength==del2use);
%         condtrix = find(cond);
%         touse = ismember(deltrix,condtrix);                     % Find which trials that are of the desired delay length are also in the desired condition
%         ix = deltrix(touse);                                    % Get the trial numbers of trials which satisfy the condition above
%         temp = squeeze(kin(sessix).dat(:,ix,featix));
%         temp = abs(temp);
%         condtemp = mySmooth(mean(temp,2,'omitnan'),41);
%         tempjaw = [tempjaw, condtemp];
%     end
% 
%     sel = tempjaw(:,1)-tempjaw(:,2);
%     delSelectivity = mean(sel(startix:stopix),'omitnan');
%     if delSelectivity<0
%         sel = -1*sel;
%     end
%     obj(sessix).jawSelectivity = mySmooth(sel,sm);
%     obj(sessix).jawvel = mean(tempjaw,2);
% end
% %% Concatenate avg jaw selectivities across all sessions
% smooth = 41;
% temp1 = [];
% blah1 = [];
% 
% delay = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent));
% go = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent));
% startix = find(obj(1).time>delay,1,'first');
% stopix = find(obj(1).time<go,1,'last');
% for sessix = 1:length(meta)
%     tempjaw = obj(sessix).jawvel;
%     temp1 = [temp1,tempjaw];
% %     temp1 = [temp1,mySmooth(tempjaw,smooth)];
%     
% end
% jv.ctrl = temp1;
% %%
% delay = 0;
% go = del2use;
% startix = find(obj(1).time>delay,1,'first');
% stopix = find(obj(1).time<go,1,'last');
% 
% temp2 = []; blah2 = [];
% for sessix = 1:length(meta)
% %     temp2 = [temp2,mySmooth(obj(sessix).jawSelectivity,smooth)];
%     temp2 = [temp2,obj(sessix).jawvel];
%     
% end
% jv.haz = temp2;
% 
% clear temp1; clear temp2
% %%
% delay = 0;
% start = params(1).tmin;
% startix = find(obj(1).time>start,1,'first');
% stopix = find(obj(1).time<delay,1,'last');
% 
% baseME = mean(jv.ctrl(startix:stopix,:),1,'omitnan');
% baseME = mean(baseME,'omitnan');
% jv.ctrl = jv.ctrl-baseME;
% 
% startix = find(obj(1).time>start,1,'first');
% stopix = find(obj(1).time<delay,1,'last');
% baseME = mean(jv.haz(startix:stopix,:),1,'omitnan');
% baseME = mean(baseME,'omitnan');
% jv.haz = jv.haz-baseME+5;     %% Hard-coded value rn to center at zero
% %% Sanity check
% % subplot(1,2,1)
% % imagesc(obj(1).time,1:16,jv.ctrl')
% % xlim([-1.3 2])
% % colorbar
% % 
% % subplot(1,2,2)
% % imagesc(obj(1).time+0.5,1:12,jv.haz')
% % xlim([-1.3 2])
% % colorbar
%  %% Sanity check
% % for i = 1:12
% %     samp = mode(obj(i).bp.ev.sample)-mode(obj(i).bp.ev.delay);
% %     trialstart = mode(obj(i).bp.ev.bitStart)-mode(obj(i).bp.ev.delay);
% %     disp(num2str(trialstart))
% % plot(obj(1).time+0.5,jv.haz(:,i)); hold on;
% % xline(0,'k--')
% % %xline(0.9,'k--')
% % xline(1.2,'k--')
% % xline(-1.3,'k--')
% % hold off;
% % pause
% % end
% %% Plot
% clearvars -except ctrlkin meta obj params jv kin meta obj params feature
% 
% colors = getColors();
% ctrlcol = colors.afc;
% hazcol = [0.5 0.5 0.5];
% go = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent));
% ctrlstop = find(obj(1).time<go,1,'last');
% del2use = 1.2;
% smooth = 10;
% alph = 0.2;
% 
% figure();
% ax = gca;
% toplot = mean(mySmooth(jv.haz,smooth),2);
% nSess = size(jv.haz,2);
% %err = 1.96*(mySmooth(std(jv.haz,0,2),smooth)./sqrt(nSess));
% err = std(mySmooth(jv.haz,smooth),0,2)./sqrt(nSess);
% shadedErrorBar(obj(1).time+0.5,toplot,err,{'Color',hazcol,'LineWidth',2},alph,ax);
% hold on;
% 
% toplot = mean(mySmooth(jv.ctrl,smooth),2);
% nSess = size(jv.ctrl,2);
% % err = 1.96*(mySmooth(std(jv.ctrl,0,2),smooth)./sqrt(nSess));
% err = std(mySmooth(jv.ctrl,smooth),0,2)./sqrt(nSess);
% shadedErrorBar(obj(1).time(1:ctrlstop),toplot(1:ctrlstop),err(1:ctrlstop),{'Color',ctrlcol,'LineWidth',2},alph,ax);
% xline(0,'LineStyle','--','Color','black','LineWidth',1.5)
% xline(go,'LineStyle','--','Color',ctrlcol,'LineWidth',1.5)
% xline(del2use,'LineStyle','--','Color',hazcol,'LineWidth',1.5)
% legend(' ',' ',' ','Haz',' ',' ',' ','Static','Delay','Go cue, static','Go cue, haz','Location','best')
% xlim([-1.3 del2use])
% xlabel('Time from delay onset (s)')
% ylabel(['Average ' feature])
% %% FUNCTIONS
% function [params,obj] = cleanUpData(params, meta,objs)
% % Reorganize params into a struct to match structure of updated data
% temp = params;
% clear params
% for sessix = 1:numel(meta)
%     temp2 = temp;
%     temp2.trialid = meta(sessix).trialid;
%     params(sessix) = temp2;
% end
% clear temp temp2
% 
% % Reorganize objs into a struct to match structure of updated data
% temp = objs;
% clear objs
% for sessix = 1:numel(meta)
%     temp2 = temp{sessix};
%     temp2.time = params(sessix).taxis; %+0.5
%     obj(sessix) = temp2;
% end
% end