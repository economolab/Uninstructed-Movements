% DECODING CDlate FROM ALL KINEMATIC FEATURES (with ridge regression,
% regularization; train/test split, cross-validation)
clear,clc,close all

whichcomp = 'LabPC';                                                % LabPC or Laptop

% Base path for code depending on laptop or lab PC
if strcmp(whichcomp,'LabPC')
    basepth = 'C:\Code';
elseif strcmp(whichcomp,'Laptop')
    basepth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code';
end

% add paths
utilspth = [basepth '\Munib Uninstruct Move\uninstructedMovements_v2'];
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
figpth = [basepth  '\Uninstructed-Movements\Fig 2'];
addpath(genpath(fullfile(figpth,'funcs')));
%% PARAMETERS
params.alignEvent          = 'firstLick'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off
params.condition(end+1) = {'R&no&~stim.enable&~autowater&~early'};              % no right, no stim, aw off
params.condition(end+1) = {'L&no&~stim.enable&~autowater&~early'};              % no left, no stim, aw off
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % all hits, no stim, aw off


params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;
params.bctype = 'reflect'; % options are : reflect  zeropad  none
%% SPECIFY DATA TO LOAD

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

meta = [];

% --- ALM ---
meta = loadJEB24_ALMVideo(meta,datapth);    

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);
%%
% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    disp(['Loading ME for session ' num2str(sessix)])
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end

% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%%
conds = [2 3];
sm = 30;

cols = getColors();

kinfeat = 'motion_energy';
kinix = find(strcmp(kin(1).featLeg,kinfeat));
for sessix = 1:length(meta)
    sesstitle = [meta(sessix).anm ' ; ' meta(sessix).date ];
    figure()
    for cc = 1:length(conds)
        switch cc
            case 1
                color = cols.rhit;
            case 2
                color = cols.lhit;

        end

        currcond = conds(cc);
        condtrix = params(sessix).trialid{currcond};
        currkin = kin(sessix).dat(:,condtrix,kinix);
        kin2plot = mySmooth(currkin,sm);
        kin2plot = mean(kin2plot,2,'omitnan');
        plot(obj(sessix).time,kin2plot,'Color',color,'LineWidth',2)
        hold on
    end
end
legend({'R ctrl','L ctrl'})
xlim([-2.5, 2.5])
xlabel(['Time from ' params(sessix).alignEvent ' (s)'])
ylabel(kinfeat)
title(sesstitle)
%% Check DLC and Bpod alignment
trix2plot = 20;
cond2plot = [2 3];

kinfeat = 'tongue_length';
featix = find(strcmp(kin(1).featLeg,kinfeat));

offset = 5;
for c = 1:length(cond2plot)
    temptrix = params(sessix).trialid{cond2plot(c)};
    subplot(1,2,c)
%     condtrix = randsample(temptrix,trix2plot);
    condtrix = temptrix;
    for i = 1:length(condtrix)
        currtrial = condtrix(i);
        if contains(kinfeat,'tongue')
            kindat = kin(sessix).dat_std(:,currtrial,featix);
        else
            kindat = mySmooth(kin(sessix).dat_std(:,currtrial,featix),smooth);
        end
        toplot = offset*i+kindat;
        plottime = obj(sessix).time;
        plot(plottime,toplot,'Color','black','LineWidth',1.1); hold on

%         allLicks = [obj(sessix).bp.ev.lickL{currtrial},obj(sessix).bp.ev.lickR{currtrial}];
%         if ~isempty(allLicks)
%             allLicks = allLicks-obj(sessix).bp.ev.goCue(currtrial);
%             plot(allLicks, max(toplot), '*', 'Color','black', 'MarkerSize',2.5);
%         end 


    rcol = cols.rhit;
    lcol = cols.lhit;

        RLicks = obj(sessix).bp.ev.lickR{currtrial};
        if ~isempty(RLicks)
            RLicks = RLicks-obj(sessix).bp.ev.(params(1).alignEvent)(currtrial);
            plot(RLicks, max(toplot), '*', 'Color',rcol, 'MarkerSize',2.5);
        end 

        LLicks = obj(sessix).bp.ev.lickL{currtrial};
        if ~isempty(LLicks)
            LLicks = LLicks-obj(sessix).bp.ev.(params(1).alignEvent)(currtrial);
            plot(LLicks, max(toplot), '*', 'Color',lcol, 'MarkerSize',2.5);
        end 
    end
    
    
    xline(0,'k--','LineWidth',1)
    title(params(sessix).condition{cond2plot(c)})
    xlim([-0.9 2.5])
    ylabel(kinfeat)

end
%%
kinfeat = 'tongue_angle';
featix = find(strcmp(kin(1).featLeg,kinfeat));

kin2plot = kin(1).dat(:,:,featix);
imagesc(obj(1).time,1:size(kin2plot,2),kin2plot')
colorbar
% clim([-0.4 0.2])