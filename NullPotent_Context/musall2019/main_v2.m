clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
% remove paths we don't need
rmpath(genpath(fullfile(utilspth,'fig1/')));
rmpath(genpath(fullfile(utilspth,'mc_stim/')));
rmpath(genpath(fullfile(utilspth,'fig3/')));

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
params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
params.condition(end+1) = {'(hit|miss)&~stim.enable&~autowater'};

params.tmin = -1.3;
params.tmax = 1.0;
params.dt = 1/200;
params.viddt = 0.0025;

% smooth with causal gaussian kernel
params.smooth = 31;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% svd params
params.reload_Cam0 = true; % redo video svd for Cam0
params.reload_Cam1 = true; % redo video svd for Cam1
params.nrDims = 200; % number of SVs to keep, 200 is probably more than needed, but that's ok

% regression params
params.frameRate = 1/params.dt / 10;
params.sPostTime = ceil(6 * params.frameRate);   % follow stim events for sPostStim in frames (used for eventType 2)
params.mPreTime = ceil(0.5 * params.frameRate);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
params.mPostTime = ceil(1 * params.frameRate);   % follow motor events for mPostStim in frames (used for eventType 3)



%% SPECIFY DATA TO LOAD

% --------------------------------------------------------- %
% TODO: make this code work for multiple sessions
% this code only works for 1 session at the moment
% --------------------------------------------------------- %
datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM ---
% meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth);
% meta = loadEKH3_ALMVideo(meta,datapth);
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);
% meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

% only keep JEB7 4-29
dates = {meta.date}';
[~,use] = patternMatchCellArray(dates,{'04-29'},'all');
meta = meta(use);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);

params.framesPerTrial = numel(obj.time);       % nr. of frames per trial


%% trials

% here, we specify trials we want to use
% this is important for the analysis, but also just for computational
% reasons
% the more trials, the more video data, and we might not have enough memory
% to load all the video in memory and perform svd
% We only load/analyze the time points we need from the video data, and we crop
% the pixels, so it *might* be possible to get away with loading all
% trials, but I haven't tested the limits yet.

params.trials2use = 1:10; % a handful of trials for testing
% params.trials2use = params.trialid{3}; % right and left, hit,miss, ~stim,~autowater
trials.right = find(ismember(params.trials2use,params.trialid{1}));
trials.left = find(ismember(params.trials2use,params.trialid{2}));



% sample equal number of left and right trials
params.nTrials2use = min(numel(trials.right,trials.left));
params.trials2use = [randsample(trials.right,params.nTrials2use,false) ; randsample(trials.left,params.nTrials2use,false)];
trials.right = sort(params.trials2use(1:params.nTrials2use));
trials.rightix = 1:params.nTrials2use;
trials.left = sort(params.trials2use((params.nTrials2use+1):end));
trials.leftix = (params.nTrials2use+1):numel(params.trials2use);


%%
clearvars -except meta params obj datapth trials

%% load raw video data if needed (or if reload set to true)

%--- Cam0 ---%
cam = 'Cam0';
[vids,vidpth,vidtm,nFrames] = getVideoMeta(cam,datapth,meta,params);

% and process if needed

% crop pixels (helps with computation and also need to remove artifacts from droplets and empty space)
xcrop = 80; % just specifying how much to crop from left end of video, will always keep from xcrop:(xwidth)
ycrop = 1:140;
[grayframes_Cam0, vidSize_Cam0] = getVideoFrames(cam,vidpth,datapth,meta,params,vids,nFrames,obj,xcrop,ycrop);


%--- Cam1 ---%
cam = 'Cam1';
[vids,vidpth,vidtm,nFrames] = getVideoMeta(cam,datapth,meta,params);

% and process if needed

% crop pixels (helps with computation and also need to remove artifacts from droplets and empty space)
xcrop = 246:371; % crop pixels (helps with computation and also need to remove artifacts from droplets and empty space)
ycrop = 65:184;
[grayframes_Cam1, vidSize_Cam1] = getVideoFrames(cam,vidpth,datapth,meta,params,vids,nFrames,obj,xcrop,ycrop);


%% check for video SVD and create if missing (or if reload set to true)

%--- Cam0 ---%
cam = 'Cam0';
sav = 0; % if 1 - save video SVs, 0 - don't save (saves to same directory as Cam0/Cam1)
vidR_Cam0 = performVideoSVD(cam,datapth,meta,params,grayframes_Cam0,vidSize_Cam0,sav); % (time,trials,singular vectors)

%--- Cam0 ---%
cam = 'Cam1';
sav = 0; % if 1 - save video SVs, 0 - don't save (saves to same directory as Cam0/Cam1)
vidR_Cam1 = performVideoSVD(cam,datapth,meta,params,grayframes_Cam1,vidSize_Cam1,sav); % (time,trials,singular vectors)


%% synchronize neural and video data

clutm = obj.time;

vidtm = params.tmin:params.viddt:params.tmax;

vidR_Cam0 = syncTimes(vidtm,vidR_Cam0,clutm);
vidR_Cam1 = syncTimes(vidtm,vidR_Cam1,clutm);


%% ready to start doing shit!!

% vidR_Cam0 are the singular vectors for cam0
% vidR_Cam1 are the singular vectors for cam0

% and then you have params, obj to get everything else you need
% and everything but obj.traj should be synced in time already


%% Helper functions

function [vids,vidpth,vidtm,nFrames] = getVideoMeta(cam,datapth,meta,params)
vidpth = fullfile(datapth,'Video',meta.anm,meta.date,cam); % assumes that data is stored in datapth/Video/<anm>/<date>/Cam0/Cam1
vids = dir(vidpth);
vids = {vids.name}';
vids = patternMatchCellArray(vids,{meta.anm,meta.date},'all');
vids = natsortfiles(vids);

vids = vids(params.trials2use); % trim trials
vidtm = params.tmin:params.viddt:params.tmax;
nFrames = numel(vidtm);
end




function [grayframes_Cam, vidSize_Cam] = getVideoFrames(cam,vidpth,datapth,meta,params,vids,nFrames,obj,xcrop,ycrop)


reload = params.(['reload_' cam]); % if true, loads video data even if SVs already exist

fpth = fullfile(datapth,'Video',meta.anm,meta.date);
if ~exist(fullfile(fpth ,['vidR_' cam   '.mat']), 'file') || ~exist(fullfile(fpth ,['absVidR_' cam   '.mat']), 'file') || reload

    for vidix = 1:numel(vids) % vidix is same as trial since vids are sorted by timestamp
        disp(['Loading vid for trial: ' num2str(vidix) '/' num2str(numel(vids))])
        vid = vids{vidix};

        % video time
        tm = obj.traj{1}(params.trials2use(vidix)).frameTimes - 0.5 - obj.bp.ev.(params.alignEvent)(params.trials2use(vidix));
        % find idx of tmin and tmax in video time
        ix1 = find(tm >= params.tmin, 1, 'first');
        ix2 = find(tm <= params.tmax, 1, 'last');
        nFramesTrial = ix2 - ix1;
        %         t1 = tm(ix1);
        %         t2 = tm(ix2);
        if nFramesTrial ~= nFrames
            ix2 = ix2 + (nFrames-nFramesTrial);
        end

        Cnt = 0;
        v = VideoReader(fullfile(vidpth, vid));
        %     cVideo = zeros(v.Width,v.Height,nFrames, 'uint8');

        frames = read(v,[ix1 ix2]);
        xcrop = xcrop:size(frames,2); % crop pixels (helps with computation and also need to remove artifacts from droplets and empty space)
        for fix = 1:nFrames
            mov = im2gray(frames(ycrop,xcrop,:,fix));
            grayframes(:,:,fix,vidix) = mov; % (x,y,frames,trials)
        end

    end
    vidSize_Cam = size(grayframes);
    grayframes_Cam = reshape(grayframes,vidSize_Cam(1),vidSize_Cam(2),vidSize_Cam(3)*vidSize_Cam(4)); % (x,y,frames*trials)
else
    disp('video and motion energy SVDs already exist and params.reload set to not reload video data - skipping')
    grayframes_Cam = nan;
    vidSize_Cam = nan;
    return
end

end


function vidR = performVideoSVD(cam,datapth,meta,params,grayframes,vidSize,sav)

reload = params.(['reload_' cam]); % if true, computes video SVs even if they already exist

fpth = fullfile(datapth,'Video',meta.anm,meta.date);

nrDims = params.nrDims; % svd dims to keep
if ~exist(fullfile(fpth, ['vidR_' cam   '.mat']), 'file') || reload
    cVideo = reshape(grayframes, [], size(grayframes,3)); %merge all pixels

    [~,s,Vr] = fsvd(single(cVideo),nrDims,1,0,0);
    vidR = s*Vr';
    vidR = bsxfun(@minus, vidR, mean(vidR, 2));
    vidR = bsxfun(@rdivide, vidR, std(vidR, [], 2))'; %this is the video regressor in the model
    vidR = reshape(vidR,vidSize(3),vidSize(4),nrDims); % (time,trials,videodims)
    if sav
        save(fullfile(fpth, ['vidR_' cam   '.mat']),'vidR');
    end
else
    temp = load(fullfile(fpth ,['vidR_' cam   '.mat']),'vidR');
    vidR = temp.vidR; clear temp;
end

end






























