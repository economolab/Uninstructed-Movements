% DECODING CDlate FROM ALL KINEMATIC FEATURES
clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Decoding Analysis'));
%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % R 2AFC hits, no stim
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % L 2AFC hits, no stim
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % R error 2AFC, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % L error 2AFC, no stim

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = (1/100)*3;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;
%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

[obj,params] = loadSessionData(meta,params);
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
