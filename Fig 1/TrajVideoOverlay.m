% DECODING CDlate FROM ALL KINEMATIC FEATURES
clear,clc,close all

whichcomp = 'Laptop';                                                % LabPC or Laptop

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
sessix = 12;
figure();
cnt = 1;
% 70:76
% 107:113
for t = 107:113
    trialnum = t;
    view = 1;
    frameTimes = obj(sessix).traj{view}(trialnum).frameTimes;
    Frames2use = find(frameTimes<(obj(sessix).bp.ev.goCue(trialnum)+2));            % Use all frames from trial start to 1s after go-cue
    allfeats = obj(sessix).traj{view}(trialnum).ts(Frames2use,2,:);
    featnames = obj(sessix).traj{view}(trialnum).featNames;
    feats2plot = [4,1];                         % Side view: 4 = jaw, 1 = tongue
    taxis = obj(sessix).traj{view}(trialnum).frameTimes(Frames2use)-obj(sessix).bp.ev.goCue(trialnum)-0.5;
    
    % Get times when a lickport contact is made
    licks = sort([obj(sessix).bp.ev.lickL{trialnum},obj(sessix).bp.ev.lickR{trialnum}],'ascend');
    licks2use = licks(find(licks<(obj(sessix).bp.ev.goCue(trialnum)+2)));
    licks2use = licks2use-obj(sessix).bp.ev.goCue(trialnum);

    % Get the times when the tongue is visible.  Plot the jaw in red when
    % the tongue is visible 
    TongueTs = allfeats(:,:,2);
    TongueVisIx = find(isnan(TongueTs));
    Jaw = allfeats(:,:,4);
    Jaw(TongueVisIx) = NaN;

    offset = 50;
    if ~obj(sessix).bp.early(trialnum)&&~obj(sessix).bp.autowater(trialnum)
        col = [0.55 0.55 0.55];
        plot(taxis,offset*cnt+(allfeats(:,:,4)),'Color',col,'LineWidth',2);
        hold on;

        col = [1 0.1 0.1];
        plot(taxis,offset*cnt+Jaw,'Color',col,'LineWidth',2);
        for l = 1:length(licks2use)
            scatter(licks2use(l),max(offset*cnt+Jaw)+5,10,'MarkerEdgeColor','black','MarkerFaceColor','black','Marker','|','LineWidth',1.75); hold on;
        end
    end
    cnt=cnt+1;
end
xline(0,'LineWidth',1.2,'LineStyle','--','Color','black')
xlim([-1 1])
sgtitle([meta(sessix).anm ' ' meta(sessix).date])
%% Plot side and bottom views of feature trajectories for pre-goCue frames
% sessix 12; trials 50:100 JEB15 07-27
sessix = 12;

figure();
for t = 50:100
    trialnum = t;
    for v = 1:2
        view = v;

        frameTimes = obj(sessix).traj{view}(trialnum).frameTimes;
        preGoFrames = find(frameTimes<obj(sessix).bp.ev.goCue(trialnum));
        allfeats = obj(sessix).traj{view}(trialnum).ts(preGoFrames,1:2,:);
        featnames = obj(sessix).traj{view}(trialnum).featNames;
        if view ==1
            feats2plot = [4,6];             % Side view: 4 = jaw, 7 = nose; Bottom view
        else
            feats2plot = [6,9,10];      % Bottom view: 5 and 6 = paws, 8 = jaw, 9 and 10 = nostrils
        end

        if ~obj(sessix).bp.early(trialnum)&&~obj(sessix).bp.autowater(trialnum)
            subplot(1,2,v)
            for f = 1:length(feats2plot)
                feat = feats2plot(f);
                if view==1
                    col = [0 0 1];
                else
                    colors = {[0 0 1], [1 0 0],[0 1 0],[1 0 1],[0 1 1]};
                    col = colors{f};
                end
                plot(allfeats(:,1,feat),allfeats(:,2,feat),'Color',col,'LineWidth',1); hold on;
            end
            set(gca, 'YDir','reverse')
            %set(gca, 'XDir','reverse')
            sgtitle([meta(sessix).anm ' ' meta(sessix).date ' ' 'Trial # ' num2str(trialnum)])
        end
    end
end
