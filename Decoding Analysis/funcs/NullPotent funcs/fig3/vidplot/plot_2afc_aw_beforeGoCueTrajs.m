% can run this script after running st_elsayed/pca/kaufman.m



clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')));
rmpath(genpath(fullfile(utilspth,'mc_stim/')));

% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'~stim.enable&~autowater&~early'};                   % afc               (1)
params.condition(end+1) = {'~stim.enable&autowater'};                    % aw                (2)
params.condition(end+1) = {'R&~stim.enable&~autowater'};             % right hit 2afc    (3)
params.condition(end+1) = {'L&~stim.enable&~autowater'};             % left hit 2afc     (4)
params.condition(end+1) = {'R&~stim.enable&autowater'};              % right hit aw      (5)
params.condition(end+1) = {'L&~stim.enable&autowater'};              % left hit aw       (6)

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

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

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM ---
% meta = loadJEB6_ALMVideo(meta,datapth);
% meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth);
% meta = loadEKH3_ALMVideo(meta,datapth);
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);


% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);



% set session data to use
anm = 'JEB15';
date = '2022-07-27';

anms = {meta(:).anm};
dates = {meta(:).date};

anmix = ismember(anms,anm);
dateix = ismember(dates,date);
sessix = find(all([anmix;dateix],1));

meta = meta(sessix);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);

% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end


%%


view = 1;

if view == 1
    vidfns = 'JEB7_2021-04-29_cam_0_date_2021_04_29_time_16_11_57_v001.avi';
else
    vidfns = 'JEB7_2021-04-29_cam_1_date_2021_04_29_time_16_11_57_v001.avi';
end


vidpth = 'C:\Users\munib\Documents\Economo-Lab\data\Video\JEB7\2021-04-29';
vidpth = fullfile(vidpth,['Cam' num2str(view-1)]);

clear v hFig hAx hIm traj
v = VideoReader(fullfile(vidpth,vidfns));

%
close all


if view == 1
    feats = [4 6];
else
    feats = [5 6 8 9 10];
end


iframe = 250;
frames = read(v,iframe);

% get obj trajectories
cond2use = [1:6];
for cix = 1:numel(cond2use)
    clear traj
    trials = randsample(params.trialid{cond2use(cix)},30,false);
    for trix = 1:numel(trials)
        trialnum = trials(trix);
        frametimes = obj.traj{view}(trialnum).frameTimes - 0.5;
        sample = obj.bp.ev.sample(trialnum);
        delay = obj.bp.ev.delay(trialnum);
        gocue = obj.bp.ev.goCue(trialnum);
        % resmaple obj.time to frametimes, but first need to get frametimes from
        % -2.5 from go cue to 2.5 to go cue
        tmin = params.tmin;
        tmax = -0.02;

%         [~,ix1] = min(abs(frametimes - gocue - tmin));
%         [~,ix2] = min(abs(frametimes - gocue - tmax));
        ix1 = 1;
        ix2 = 980; % just before go cue
%         ix2 - ix1

        frametimes = frametimes(ix1:ix2);

        % trajectories across all trials, before go cue

        
        traj(:,:,:,trix) = obj.traj{view}(trialnum).ts(ix1:ix2,[1 2],feats);
        
    end
    traj = permute(traj, [1 4 2 3]);
    traj = reshape(traj,size(traj,1)*size(traj,2),size(traj,3),size(traj,4));

    cols = linspecer(numel(feats),'qualitative');

    f = figure;
    hold on;
    imagesc(flipud(frames));

    for fix = 1:numel(feats)
        xs = round(squeeze(traj(:,1,fix)));
        ys = v.Height - round(squeeze(traj(:,2,fix)));

        ff = plot(xs, ys , 'Color', cols(fix,:),'LineWidth',1);
        ff.Color = [ff.Color 0.05];
    end
    xlim([0 v.Width])
    ylim([0 v.Height])
    title(params.condition{cond2use(cix)},'Interpreter','none')

end
















