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
%% data object and video file
figure();
sessix = 3;
trialrange = 50:100;
for view = 1:2                                                                         % For side cam and bottom cam...
    if view==1
        datapth = 'C:\Users\Jackie Birnbaum\Documents\Data\DataObjects\JEB13';         % Directory where data obj is stored
        objfn = 'data_structure_JEB13_2022-09-25';                                     % File name of data object
        vidfn = 'JEB13_2022-09-25_cam_0_date_2022_09_25_time_09_52_35_v001.avi';       % Filename of video
        v = VideoReader(fullfile(datapth, vidfn));
    else
        datapth = 'C:\Users\Jackie Birnbaum\Documents\Data\DataObjects\JEB13';         % Directory where data obj is stored
        objfn = 'data_structure_JEB13_2022-09-25';                                     % File name of data object
        vidfn = 'JEB13_2022-09-25_cam_1_date_2022_09_25_time_09_52_35_v001.avi';       % Filename of video
        v = VideoReader(fullfile(datapth, vidfn));
    end
    % Plot side and bottom views of feature trajectories for pre-goCue frames
    subplot(1,2,view)
    imagesc(readFrame(v)); hold on;                                                    % Plot camera image

    probThresh = 0.9999;                                                               % Remove trajectory points that have a DLC confidence level lower than this threshold 
    for t = trialrange                                                                 % Trial range
        trialnum = t;
        frameTimes = obj(sessix).traj{view}(trialnum).frameTimes;                      % Get the frame times for this trial
        preGoFrames = find(frameTimes<obj(sessix).bp.ev.goCue(trialnum));              % Take frames that come before the goCue
        allfeats = obj(sessix).traj{view}(trialnum).ts(preGoFrames,:,:);               % Get all of the feature trajectories for pre-goCue frames
        featnames = obj(sessix).traj{view}(trialnum).featNames;                        % Feature names
        if view ==1
            feats2plot = [4,6];             % Side view: 4 = jaw, 6 = nose             % Which features you want to plot from appropriate camera
        else
            feats2plot = [5,6,8,9,10];          % Bottom view: 5 and 6 = paws, 8 = jaw, 9 and 10 = nostrils
        end

        if ~obj(sessix).bp.early(trialnum)&&~obj(sessix).bp.autowater(trialnum)        % If this trial is not an autowater and not an early lick trial  
            for f = 1:length(feats2plot)                                               % For each feature
                feat = feats2plot(f);
                if view==1
                    col = [0 0 1];
                else
                    colors = {[0 0 1], [1 0 0],[0 1 0],[1 0 1],[0 1 1]};
                    col = colors{f};
                end
                outlierFrames = find(allfeats(:,3,feat)<probThresh);                   % Get the indices of frames where the current feature has a confidence level below the threshold
                if ~isempty(outlierFrames)                                             % If there are any outlier frames...
                    allfeats(outlierFrames,:,feat) = NaN;                              % Get rid of them
                end
                plot(allfeats(:,1,feat),allfeats(:,2,feat),'Color',col,'LineWidth',1); hold on;     % Plot feature trajectories overlaid onto a still frame from video
            end
            set(gca, 'YDir','reverse')
            %set(gca, 'XDir','reverse')
            sgtitle([meta(sessix).anm ' ' meta(sessix).date ' ' 'Trial # ' num2str(trialnum)])
        end
    end
end