% DECODING CDlate FROM ALL KINEMATIC FEATURES
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
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'~autowater&~early&~no'};                         % R 2AFC hits, no stim
params.condition(end+1) = {'autowater&~early&~no'};                         % R 2AFC hits, no stim
params.condition(end+1) = {'~early&~no'};                         % R 2AFC hits, no stim

params.tmin = -3;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril','top_paw','bottom_paw'}};

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
meta = loadJEB13_ALMVideo(meta,datapth);
% meta = loadJEB6_ALMVideo(meta,datapth);
% meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth);
% meta = loadEKH3_ALMVideo(meta,datapth);
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);
% meta = loadJEB15_ALMVideo(meta,datapth);

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
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%%
times.start = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent));
times.startix = find(obj(1).time>times.start,1,'first');
times.stop = mode(obj(1).bp.ev.sample)-mode(obj(1).bp.ev.(params(1).alignEvent))-0.05;
times.stopix = find(obj(1).time<times.stop,1,'last');

cond2use = 4;

feat2use = {'jaw_yvel_view1','nose_yvel_view1', 'top_paw_yvel_view2'};
featix = NaN(1,length(feat2use));
for f = 1:length(feat2use)
    currfeat = feat2use{f};
    currix = find(strcmp(kin(1).featLeg,currfeat));
    featix(f) = currix;
end

trix2use = 100;
% load([basepth '\Uninstructed-Movements\Fig 1\RedColormap.mat']);
% load([basepth '\Uninstructed-Movements\Fig 1\BlueColormap.mat']);
% load([basepth '\Uninstructed-Movements\Fig 1\GreenColormap.mat']);
% cmaps = {RedColormap,BlueColormap,GreenColormap};


del = median(obj(1).bp.ev.delay)-median(obj(1).bp.ev.(params(1).alignEvent));
delix = find(obj(1).time>del,1,'first');
go = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent));
goix = find(obj(1).time<go,1,'last');

sm = 61;


for sessix = 5
    figure();
    allkin = [];
    cond = cond2use;
    condtrix = params(sessix).trialid{cond};
    ntrials = length(condtrix);
    randtrix = randsample(condtrix,trix2use);
    for f = 1:length(featix)
        currfeat = featix(f);
        if f~=3
            currkin = mySmooth(kin(sessix).dat_std(times.startix:goix,randtrix,currfeat),sm);
            currkin = abs(currkin);
        else
            % Mean-center paw
            tempkin = kin(sessix).dat_std(:,randtrix,currfeat);
            presampME = squeeze(mean(tempkin(times.startix:times.stopix,:,:),1,'omitnan'));
            avgpresampME = mean(presampME,'omitnan');
            tempkin = tempkin-avgpresampME;
            currkin = mySmooth(tempkin,sm);
            currkin = currkin(times.startix:goix,:);
        end

        % Max normalize current feature
        maxkin = max(currkin,[],"all");
        currkin = abs(currkin./maxkin);

        

        allkin = cat(3,allkin,currkin);
    end

    allkin = permute(allkin,[2 1 3]);
    IM = imshow(allkin,'InitialMagnification','fit');
    title(['RGB = ' feat2use])
end
%%
onemat = ones(size(allkin));
newkin = onemat-allkin;

RI = imref2d(size(newkin));
RI.XWorldLimits = [0 3];
RI.YWorldLimits = [2 5];
IMref = imshow(newkin, RI,'InitialMagnification','fit');
title(['RGB = ' feat2use])
% for tt = 1:trix2use
%     currtrial = squeeze(allkin(:,tt,:));
%     for tw = 1:size(allkin,1)
%         currwin = currtrial(tw,:);
% 
%         disp('hi')
%     end
% end
%%
RI = imref2d(size(allkin));
RI.XWorldLimits = [0 3];
RI.YWorldLimits = [2 5];
IMref = imshow(allkin, RI,'InitialMagnification','fit');
title(['RGB = ' feat2use])