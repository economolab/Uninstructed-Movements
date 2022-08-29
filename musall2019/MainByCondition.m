clear,clc,close all
addpath(genpath(pwd))
addpath(genpath('C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\DataLoadingScripts'));
addpath(genpath('C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\funcs'));
addpath(genpath('C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\utils'));

%% SET RUN PARAMS
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
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'all'};


% regression params
params.frameRate = 1/params.dt / 10;
params.sPostTime = ceil(6 * params.frameRate);   % follow stim events for sPostStim in frames (used for eventType 2)
params.mPreTime = ceil(0.5 * params.frameRate);  % precede motor events to capture preparatory activity in frames (used for eventType 3)
params.mPostTime = ceil(1 * params.frameRate);   % follow motor events for mPostStim in frames (used for eventType 3)
params.reload_Cam0 = false; % redo video svd
params.reload_Cam1 = false; % redo video svd

%% SET METADATA

datapth = 'C:\Users\munib\Documents\Economo-Lab\data';

meta = [];
meta = loadJEB7_ALMVideo(meta,datapth); % done

% only keep 4-29
dates = {meta.date}';
[~,use] = patternMatchCellArray(dates,{'04-29'},'all');
meta = meta(use);

params.probe = [meta.probe]; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

%% LOAD AND PROCESS DATA

objs = loadObjs(meta);


for metaix = 1:numel(meta)
    obj = objs{metaix};
    disp('______________________________________________________')
    disp(['Processing data for session ' [meta(metaix).anm '_' meta(metaix).date]])
    disp(' ')
    [sessparams{metaix},sessobj{metaix}] = processData(obj,params,params.probe(metaix));
end

% clean up sessparams and sessobj
for metaix = 1:numel(meta)
    params(metaix).trialid = sessparams{metaix}.trialid;
    params(metaix).cluid{metaix} = sessparams{metaix}.cluid{params.probe(metaix)};

    objs{metaix} = sessobj{metaix};
    objs{metaix}.psth = objs{metaix}.psth{params.probe(metaix)};
    objs{metaix}.trialdat = objs{metaix}.trialdat{params.probe(metaix)};
    objs{metaix}.presampleFR = objs{metaix}.presampleFR{params.probe(metaix)};
    objs{metaix}.presampleSigma = objs{metaix}.presampleSigma{params.probe(metaix)};
end

obj = objs{1};
params.framesPerTrial = numel(obj.time);       % nr. of frames per trial

% trials
params.trials2use = params.trialid{3}; % right and left, hit,miss, ~stim,~autowater
trials.right = find(ismember(params.trials2use,params.trialid{1}));
trials.left = find(ismember(params.trials2use,params.trialid{2}));
% sample equal number of left and right trials
params.nTrials2use = min(numel(trials.right,trials.left));
params.trials2use = [randsample(trials.right,params.nTrials2use,false) ; randsample(trials.left,params.nTrials2use,false)];
trials.right = sort(params.trials2use(1:params.nTrials2use));
trials.rightix = 1:params.nTrials2use;
trials.left = sort(params.trials2use((params.nTrials2use+1):end));
trials.leftix = (params.nTrials2use+1):numel(params.trials2use);

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')


%%
clearvars -except meta params obj datapth trials

%% load raw video data and process if needed

%--- Cam0 ---%
cam = 'Cam0';

vidpth = fullfile(datapth,'Video',meta.anm,meta.date,cam);
vids = dir(vidpth);
vids = {vids.name}';

vids = patternMatchCellArray(vids,{meta.anm,meta.date},'all');
vids = natsortfiles(vids);

vids = vids(params.trials2use); % trim trials

vidtm = params.tmin:params.viddt:params.tmax;
nFrames = numel(vidtm);

fpth = fullfile(datapth,'Video',meta.anm,meta.date);
if ~exist(fullfile(fpth ,['vidR_' cam   '.mat']), 'file') || ~exist(fullfile(fpth ,['absVidR_' cam   '.mat']), 'file') || params.reload_Cam0

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
        xcrop = 80:size(frames,2); % crop pixels (helps with computation and also need to remove artifacts from droplets and empty space)
        ycrop = 1:140;
        for fix = 1:nFrames
            %         grayframes(:,:,fix) = im2gray(frames(:,:,:,fix));
            mov = im2gray(frames(ycrop,xcrop,:,fix));
            grayframes(:,:,fix,vidix) = mov; % (x,y,frames,trials)
        end

    end
    vidSize_Cam0 = size(grayframes);
    grayframes_Cam0 = reshape(grayframes,vidSize_Cam0(1),vidSize_Cam0(2),vidSize_Cam0(3)*vidSize_Cam0(4)); % (x,y,frames*trials)
end

%--- Cam1 ---%
cam = 'Cam1';

clear grayframes

vidpth = fullfile(datapth,'Video',meta.anm,meta.date,cam);
vids = dir(vidpth);
vids = {vids.name}';

vids = patternMatchCellArray(vids,{meta.anm,meta.date},'all');
vids = natsortfiles(vids);

vids = vids(params.trials2use); % trim trials

vidtm = params.tmin:params.viddt:params.tmax;
nFrames = numel(vidtm);

fpth = fullfile(datapth,'Video',meta.anm,meta.date);
if ~exist(fullfile(fpth ,['vidR_' cam   '.mat']), 'file') || ~exist(fullfile(fpth ,['absVidR_' cam   '.mat']), 'file') || params.reload_Cam1

    for vidix = 1:numel(vids) % vidix is same as trial since vids are sorted by timestamp
        disp(['Loading vid for trial: ' num2str(vidix) '/' num2str(numel(vids))])
        vid = vids{vidix};

        % video time
        tm = obj.traj{2}(params.trials2use(vidix)).frameTimes - 0.5 - obj.bp.ev.(params.alignEvent)(params.trials2use(vidix));
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
        xcrop = 246:371; % crop pixels (helps with computation and also need to remove artifacts from droplets and empty space)
        ycrop = 65:184;
        for fix = 1:nFrames
            %         grayframes(:,:,fix) = im2gray(frames(:,:,:,fix));
            mov = im2gray(frames(ycrop,xcrop,:,fix));
            grayframes(:,:,fix,vidix) = mov; % (x,y,frames,trials)
        end

    end
    vidSize_Cam1 = size(grayframes);
    grayframes_Cam1 = reshape(grayframes,vidSize_Cam1(1),vidSize_Cam1(2),vidSize_Cam1(3)*vidSize_Cam1(4)); % (x,y,frames*trials)
end



%% check for video SVD and create if missing

%--- Cam0 ---%
cam = 'Cam0';
fpth = fullfile(datapth,'Video',meta.anm,meta.date);

nrDims = 200; % svd dims to keep
if ~exist(fullfile(fpth, ['vidR_' cam   '.mat']), 'file') || params.reload_Cam0
    cVideo = reshape(grayframes_Cam0, [], size(grayframes_Cam0,3)); %merge all pixels

    [~,s,Vr] = fsvd(single(cVideo),nrDims,1,0,0);
    vidR = s*Vr';
    vidR = bsxfun(@minus, vidR, mean(vidR, 2));
    vidR = bsxfun(@rdivide, vidR, std(vidR, [], 2))'; %this is the video regressor in the model
    vidR = reshape(vidR,vidSize_Cam0(3),vidSize_Cam0(4),nrDims); % (time,trials,videodims)
    save(fullfile(fpth, ['vidR_' cam   '.mat']),'vidR');
else
    temp = load(fullfile(fpth ,['vidR_' cam   '.mat']),'vidR');
    vidR_Cam0 = temp.vidR; clear temp;
end

%--- Cam1 ---%
cam = 'Cam1';
fpth = fullfile(datapth,'Video',meta.anm,meta.date);

nrDims = 200; % svd dims to keep
if ~exist(fullfile(fpth, ['vidR_' cam   '.mat']), 'file') || params.reload_Cam1
    cVideo = reshape(grayframes_Cam1, [], size(grayframes_Cam1,3)); %merge all pixels

    [U,s,Vr] = fsvd(single(cVideo),nrDims,1,0,0);
    vidR = s*Vr';
    vidR = bsxfun(@minus, vidR, mean(vidR, 2));
    vidR = bsxfun(@rdivide, vidR, std(vidR, [], 2))'; %this is the video regressor in the model
    vidR = reshape(vidR,vidSize_Cam1(3),vidSize_Cam1(4),nrDims); % (time,trials,videodims)
    save(fullfile(fpth, ['vidR_' cam   '.mat']),'vidR');
else
    temp = load(fullfile(fpth ,['vidR_' cam   '.mat']),'vidR');
    vidR_Cam1 = temp.vidR; clear temp;
end



%% compute motion energy and repeat svd if needed

%--- Cam0 ---%
cam = 'Cam0';
fpth = fullfile(datapth,'Video',meta.anm,meta.date);

if ~exist(fullfile(fpth ,['absVidR_' cam   '.mat']), 'file') || params.reload_Cam0
    cVideo = grayframes_Cam0; %merge all pixels
    cVideo = cat(3,cVideo(:,:,1),cVideo); %duplicate first frame
    cVideo = abs(diff(cVideo,[],3)); %compute absolute of temporal derivative (motion energy)
    cVideo = reshape(cVideo, [], size(cVideo,3)); %merge all pixels

    [~, s, Vr] = fsvd(single(cVideo), nrDims, 1, 0, 0);
    absVidR = s * Vr';
    absVidR = bsxfun(@minus, absVidR, mean(absVidR, 2));
    absVidR = bsxfun(@rdivide, absVidR, std(absVidR, [], 2))'; %this is the motion energy regressor in the model
    absVidR = reshape(absVidR,vidSize_Cam0(3),vidSize_Cam0(4),nrDims); % (time,trials,videodims)
    save(fullfile(fpth, ['absVidR_' cam   '.mat']),'absVidR');
else
    temp = load(fullfile(fpth ,['absVidR_' cam   '.mat']),'absVidR');
    absVidR_Cam0 = temp.absVidR; clear temp;
end
clear cVideo

%--- Cam1 ---%
cam = 'Cam1';
fpth = fullfile(datapth,'Video',meta.anm,meta.date);

if ~exist(fullfile(fpth ,['absVidR_' cam   '.mat']), 'file') || params.reload_Cam1
    cVideo = grayframes_Cam1; %merge all pixels
    cVideo = cat(3,cVideo(:,:,1),cVideo); %duplicate first frame
    cVideo = abs(diff(cVideo,[],3)); %compute absolute of temporal derivative (motion energy)
    cVideo = reshape(cVideo, [], size(cVideo,3)); %merge all pixels

    [~, s, Vr] = fsvd(single(cVideo), nrDims, 1, 0, 0);
    absVidR = s * Vr';
    absVidR = bsxfun(@minus, absVidR, mean(absVidR, 2));
    absVidR = bsxfun(@rdivide, absVidR, std(absVidR, [], 2))'; %this is the motion energy regressor in the model
    absVidR = reshape(absVidR,vidSize_Cam1(3),vidSize_Cam1(4),nrDims); % (time,trials,videodims)
    save(fullfile(fpth, ['absVidR_' cam   '.mat']),'absVidR');
else
    temp = load(fullfile(fpth ,['absVidR_' cam   '.mat']),'absVidR');
    absVidR_Cam1 = temp.absVidR; clear temp;
end
clear cVideo


%% synchronize neural and video data

clutm = obj.time;

vidtm = params.tmin:params.viddt:params.tmax;

vidR_Cam0 = syncTimes(vidtm,vidR_Cam0,clutm);
vidR_Cam1 = syncTimes(vidtm,vidR_Cam1,clutm);

absVidR_Cam0 = syncTimes(vidtm,absVidR_Cam0,clutm);
absVidR_Cam1 = syncTimes(vidtm,absVidR_Cam1,clutm);

%% kinematics for comparison with bottom cam svd

% dfparams.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
%                         {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};
% dfparams.advance_movement = 0.0;
% [kin,kinfeats] = getKin(meta,obj,dfparams,params);
% 

%% ephys data

N = permute(obj.trialdat,[1 3 2]); % single trial binned neural data (time,trials,clusters)

% trim trials
N = N(:,params.trials2use,:);

sz.N = size(N);
sz.vidR_Cam0 = size(vidR_Cam0);
sz.vidR_Cam1 = size(vidR_Cam1);
sz.absVidR_Cam0 = size(absVidR_Cam0);
sz.absVidR_Cam1 = size(absVidR_Cam1);
sz.leg = {'time','trials','clu/dims'};

%% video data
% instructed video is average post-gocue video 
% uninstructed video is the remaining (variability that is not captured
% with trial-averaged data)

[~,goCueIx] = min(abs(obj.time - 0));
toPlot = 0;

ix = 1; % right trials
[instructed_vidR_Cam0{ix}, uninstructed_vidR_Cam0{ix}]       = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.rightix,:),goCueIx,toPlot);
[instructed_vidR_Cam1{ix}, uninstructed_vidR_Cam1{ix}]       = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.rightix,:),goCueIx,toPlot);
[instructed_absVidR_Cam0{ix}, uninstructed_absVidR_Cam0{ix}] = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.rightix,:),goCueIx,toPlot);
[instructed_absVidR_Cam1{ix}, uninstructed_absVidR_Cam1{ix}] = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.rightix,:),goCueIx,toPlot);

ix = 2; % left trials
[instructed_vidR_Cam0{ix}, uninstructed_vidR_Cam0{ix}]       = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.leftix,:),goCueIx,toPlot);
[instructed_vidR_Cam1{ix}, uninstructed_vidR_Cam1{ix}]       = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.leftix,:),goCueIx,toPlot);
[instructed_absVidR_Cam0{ix}, uninstructed_absVidR_Cam0{ix}] = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.leftix,:),goCueIx,toPlot);
[instructed_absVidR_Cam1{ix}, uninstructed_absVidR_Cam1{ix}] = getInstructedAndUninstructedSVDs(vidR_Cam0(:,trials.leftix,:),goCueIx,toPlot);


%% reshape video and ephys data

instructed_vidR_Cam0 = catTimeTrialsSVDs(instructed_vidR_Cam0);
uninstructed_vidR_Cam0 = catTimeTrialsSVDs(uninstructed_vidR_Cam0);
instructed_vidR_Cam1 = catTimeTrialsSVDs(instructed_vidR_Cam1);
uninstructed_vidR_Cam1 = catTimeTrialsSVDs(uninstructed_vidR_Cam1);
instructed_absVidR_Cam0 = catTimeTrialsSVDs(instructed_absVidR_Cam0);
uninstructed_absVidR_Cam0 = catTimeTrialsSVDs(uninstructed_absVidR_Cam0);
instructed_absVidR_Cam1 = catTimeTrialsSVDs(instructed_absVidR_Cam1);
uninstructed_absVidR_Cam1 = catTimeTrialsSVDs(uninstructed_absVidR_Cam1);


% figure; imagesc((1:size(vidR,1)).*params.dt,1:size(vidR,2),vidR'); colorbar; caxis([-4 4])
% figure; imagesc((1:size(absVidR,1)).*params.dt,1:size(absVidR,2),absVidR'); colorbar; caxis([-4 4])

N = reshape(N,size(N,1)*size(N,2),size(N,3)); % (time*trials,clu)
N = bsxfun(@minus, N, mean(N)); %make zero-mean

% figure; imagesc((1:size(N,1)).*params.dt,1:size(N,2),N'); colorbar

%% task regressors

%-- choice
% if licked right = 0, left = 1
choice = (obj.bp.hit & obj.bp.L) | (obj.bp.miss & obj.bp.R);
choice = choice(params.trials2use);

%-- previous choice (every trial after a left choice trial)
% if licked right = 0, left = 1
prevchoice = [0 ; choice(1:end-1)];

%-- success
success = obj.bp.hit(params.trials2use);

%-- previous success
prevsuccess = [0 ; success(1:end-1)];

%-- water given (all frames after reward given)
reward = obj.bp.ev.reward - obj.bp.ev.goCue;
reward = reward(params.trials2use);
rewardix = nan(size(reward));
for i = 1:numel(reward)
    if ~isnan(reward(i))
        [~,rewardix(i)] = min(abs(obj.time - reward(i)));
    end
end


%-- now shape them accordingly

% choice
temp = zeros(numel(choice),numel(obj.time));
temp(:,1) = choice;
choice = temp'; % (time,trials)
choice = choice(:);
taskEvents(:,1) = choice;

% previous choice
temp = zeros(numel(prevchoice),numel(obj.time));
temp(:,1) = prevchoice;
prevchoice = temp'; % (time,trials)
prevchoice = prevchoice(:);
taskEvents(:,2) = prevchoice;

% success
temp = zeros(numel(success),numel(obj.time));
temp(:,1) = success;
success = temp'; % (time,trials)
success = success(:);
taskEvents(:,3) = success;

% previous success
temp = zeros(numel(prevsuccess),numel(obj.time));
temp(:,1) = prevsuccess;
prevsuccess = temp'; % (time,trials)
prevsuccess = prevsuccess(:);
taskEvents(:,4) = prevsuccess;

% reward
temp = zeros(numel(reward),numel(obj.time));
for i = 1:numel(rewardix)
    if ~isnan(rewardix(i))
        temp(i,rewardix(i)) = 1;
    end
end
reward = temp';
reward = reward(:);
taskEvents(:,5) = reward;

% figure; imagesc((1:size(taskEvents,1)).*params.dt,1:size(taskEvents,2),taskEvents')

%% make task design matrix



taskLabels = {'choice' 'prevchoice' 'success' 'prevsuccess' 'reward'}; %some task variables
taskEventType = [1 1 1 1 1]; %different type of events.

[taskR, taskIdx] = makeDesignMatrix(taskEvents, taskEventType, params); %make design matrix for task variables

% rewardR = taskR(:,taskIdx==5);


%% movement regressors

% just using licks. licks before go cue are uninstructed licks
% licks after go cue are instructed licks

licks.L = obj.bp.ev.lickL(params.trials2use);
licks.R = obj.bp.ev.lickR(params.trials2use);
lick = zeros(numel(params.trials2use),numel(obj.time));
for trix = 1:numel(params.trials2use)
    L = licks.L{trix} - obj.bp.ev.goCue(params.trials2use(trix));
    R = licks.R{trix} - obj.bp.ev.goCue(params.trials2use(trix));
    
    L = L(L>=params.tmin & L<=params.tmax);
    R = R(R>=params.tmin & R<=params.tmax);

    licktm = sort([L R])';
    for i = 1:numel(licktm)
        [~,ix] = min(abs(obj.time - licktm(i)));
        lick(trix,ix) = 1;
    end
end
lick = lick'; % (time, trials) binary matrix, 1 if lick at that time, 0 o/w

moveEvents(:,1) = lick(:);


%% make movement design matrix

moveLabels = {'licks'}; % some movement variables
moveEventType = [3]; %different type of events.

[moveR, moveIdx] = makeDesignMatrix(moveEvents, moveEventType, params); %make design matrix for task variables


% figure; plot(moveR(:,moveIdx==1),'r')
% hold on
% plot(moveR(:,moveIdx==2),'g')


%% make full design matrix

% not using instructed_vidR since it's rank-deficient

fullR = [taskR, moveR, ...
         uninstructed_absVidR_Cam0, instructed_absVidR_Cam0, ...
         uninstructed_absVidR_Cam1, instructed_absVidR_Cam1, ...
         uninstructed_vidR_Cam0, instructed_vidR_Cam0, ...
         uninstructed_vidR_Cam1, instructed_vidR_Cam1]; %make new, single design matrix
% fullR = [taskR, moveR, instructed_absVidR, uninstructed_absVidR, uninstructed_vidR]; %make new, single design matrix


fullR = bsxfun(@minus, fullR, mean(fullR, 1));

moveLabels = [moveLabels, {'uiME_Cam0'}, {'iME_Cam0'}, {'uiME_Cam1'}, {'iME_Cam1'},  {'uiVideo_Cam0'}, {'iVideo_Cam0'},  {'uiVideo_Cam1'}, {'iVideo_Cam1'}];
% moveLabels = [moveLabels, {'iME'}, {'uiME'}, {'uiVideo'}];


% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
regLabels = cat(2,taskLabels,moveLabels);

% % index to reconstruct different response kernels
regIdx = [
            ones(sum(taskIdx==1),1) * find(ismember(regLabels,'choice')); ...
            ones(sum(taskIdx==2),1) * find(ismember(regLabels,'prevchoice')); ...
            ones(sum(taskIdx==3),1) * find(ismember(regLabels,'success')); ...
            ones(sum(taskIdx==4),1) * find(ismember(regLabels,'prevsuccess')); ...
            ones(sum(taskIdx==5),1) * find(ismember(regLabels,'reward')); ...
            ones(sum(moveIdx==1),1) * find(ismember(regLabels,'licks')); ...
            ones(nrDims,1) * find(ismember(regLabels,'uiME_Cam0')); ...
            ones(nrDims,1) * find(ismember(regLabels,'iME_Cam0')); ...
            ones(nrDims,1) * find(ismember(regLabels,'uiME_Cam1')); ...
            ones(nrDims,1) * find(ismember(regLabels,'iME_Cam1')); ...
            ones(nrDims,1) * find(ismember(regLabels,'uiVideo_Cam0')); ...
            ones(nrDims,1) * find(ismember(regLabels,'iVideo_Cam0')); ...
            ones(nrDims,1) * find(ismember(regLabels,'uiVideo_Cam1')); ...
            ones(nrDims,1) * find(ismember(regLabels,'iVideo_Cam1')); ...
         ];

% regIdx = [
%             ones(sum(taskIdx==1),1) * find(ismember(regLabels,'choice')); ...
%             ones(sum(taskIdx==2),1) * find(ismember(regLabels,'prevchoice')); ...
%             ones(sum(taskIdx==3),1) * find(ismember(regLabels,'success')); ...
%             ones(sum(taskIdx==4),1) * find(ismember(regLabels,'prevsuccess')); ...
%             ones(sum(taskIdx==5),1) * find(ismember(regLabels,'reward')); ...
%             ones(sum(moveIdx==1),1) * find(ismember(regLabels,'licks')); ...
%             ones(nrDims*2,1) * find(ismember(regLabels,'iME')); ...
%             ones(nrDims,1) * find(ismember(regLabels,'uiME')); ...
%             ones(nrDims*2,1) * find(ismember(regLabels,'iVideo')); ...
%             ones(nrDims,1) * find(ismember(regLabels,'uiVideo')) ...
%          ];



%% run QR and check for rank-defficiency. This will show whether a given regressor is highly collinear with other regressors in the design matrix.
% The resulting plot ranges from 0 to 1 for each regressor, with 1 being
% fully orthogonal to all preceeding regressors in the matrix and 0 being
% fully redundant. Having fully redundant regressors in the matrix will
% break the model, so in this example those regressors are removed. In
% practice, you should understand where the redundancy is coming from and
% change your model design to avoid it in the first place!

rejIdx = false(1,size(fullR,2));
[~, fullQRR] = qr(bsxfun(@rdivide,fullR,sqrt(sum(fullR.^2))),0); %orthogonalize normalized design matrix
figure; plot(abs(diag(fullQRR)),'linewidth',2); ylim([0 1.1]); title('Regressor orthogonality'); drawnow; %this shows how orthogonal individual regressors are to the rest of the matrix
axis square; ylabel('Norm. vector angle'); xlabel('Regressors');
if sum(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1))) < size(fullR,2) %check if design matrix is full rank
    temp = ~(abs(diag(fullQRR)) > max(size(fullR)) * eps(fullQRR(1)));
    fprintf('Design matrix is rank-defficient. Removing %d/%d additional regressors.\n', sum(temp), sum(~rejIdx));
    rejIdx(~rejIdx) = temp; %reject regressors that cause rank-defficint matrix
end

%% run full model fit for single trials
[~, dimBeta] = ridgeMML(N, fullR, true); %make model fit
fullFit = fullR * dimBeta; %fit data

disp('Full model fit completed');


%%

uninstructedLabels = {'uiME_Cam0','uiVideo_Cam0','uiME_Cam1','uiVideo_Cam1'};
lbIdx = find(ismember(regLabels,uninstructedLabels));
uiIdx = ismember(regIdx ,lbIdx) ; % logical of uninstructed vars in cols of fullR

instructedLabels = {'licks','iME_Cam0','iVideo_Cam0','iME_Cam1','iVideo_Cam1'};
lbIdx = find(ismember(regLabels,instructedLabels));
iIdx = ismember(regIdx ,lbIdx) ; % logical of uninstructed vars in cols of fullR

lbIdx = find(ismember(regLabels,taskLabels));
taskIdx = ismember(regIdx ,lbIdx) ; % logical of uninstructed vars in cols of fullR

ui_fit = fullR(:, uiIdx) * dimBeta(uiIdx, :); % uninstructed vars fit 
i_fit = fullR(:,iIdx) * dimBeta(iIdx,:);      % instructed vars fit 
task_fit = fullR(:,taskIdx) * dimBeta(taskIdx,:);   % instructed vars fit 


%% variance explained (https://economictheoryblog.com/2014/11/05/proof/)
R = corrcoef(N,fullFit);
ve_full = R(1,2)^2;

R = corrcoef(N,i_fit);
ve_i = R(1,2)^2;

R = corrcoef(N,ui_fit);
ve_ui = R(1,2)^2;

R = corrcoef(N,task_fit);
ve_task = R(1,2)^2;

ve = [ve_full,ve_ui,ve_i,ve_task];

%
figure;
X = categorical({'Full','Unin. Move','In. Move','Task'});
h = histogram('Categories', X, 'BinCounts', ve, 'EdgeColor','none');
ylabel('Frac. VE')
ax = gca;
ax.FontSize = 13;


%%


figure; 
subplot(2,1,1)
imagesc((1:size(N,1)).*params.dt,1:size(N,2),N'); 
colorbar
caxis([-50 100])
title('Original')
xlabel('Time (s)')
ylabel('Neurons')

subplot(2,1,2)
imagesc((1:size(N,1)).*params.dt,1:size(N,2),fullFit'); 
colorbar
caxis([-50 100])
title('Reconstructed')
xlabel('Time (s)')
ylabel('Neurons')

sgtitle([num2str(ve_full*100) '% VE '])

%% example trial
figure; 
subplot(2,1,1)
imagesc(N(1000:2000,:)'); 
colorbar
caxis([-50 100])
title('Original')
xlabel('Time (s)')
ylabel('Neurons')

subplot(2,1,2)
imagesc(fullFit(1000:2000,:)'); 
colorbar
caxis([-50 100])
title('Reconstructed')
% xlabel('Time (s)')
ylabel('Neurons')

sgtitle([num2str(ve_full*100) '% VE '])


%% reconstruct PSTHs full fit

% N is original data in (time*trials,neurons)
% fullFit is reocnstructed data of same size

N = reshape(N,sz.N(1),sz.N(2),sz.N(3));
fullFit = reshape(fullFit,sz.N(1),sz.N(2),sz.N(3));

% pick a random neuron
k = 10;

%%
clus = randsample(size(N,3),k,false);
f = figure;
f.Position = [203         332        1418         635];
t = tiledlayout('flow');
for i = 1:k
    clu = clus(i);

    psthN = mean(N(:,:,clu),2);
    psthFullFit = mean(fullFit(:,:,clu),2);


    ax = nexttile; hold on;
    patchline(obj.time,psthN,'EdgeColor',[0 0 0],'LineWidth',5,'EdgeAlpha',0.3);
    patchline(obj.time,psthFullFit,'EdgeColor',[41, 128, 64]./255,'LineWidth',2,'EdgeAlpha',1)
    xlim([obj.time(10) obj.time(end)])
    ax.FontSize = 13;
end
title(t, 'PSTHs and full fit reconstructions')
xlabel(t, 'Time (s) from go cue')
ylabel(t, 'Centered Firing Rate (Hz)')

pth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\musall2019\figs';
fn = 'psth_fullfit_recon';
mysavefig(f,pth,fn)


%% reconstruct PSTHs and few single trials, task fit

% N is original data in (time*trials,neurons)
% fullFit is reocnstructed data of same size

task_fit = reshape(task_fit,sz.N(1),sz.N(2),sz.N(3));

%
clu = randsample(size(N,3),1); % sample 1 neuron
trix = sort(randsample(size(N,2),3,false)); % sample 3 trials

% clu = 19;
% trix = sort([135 138 171]);

sav = 0;
gclw = 2;
sm = 11;

close all
f1 = figure;
f1.Position = [477    82   379   897];
t = tiledlayout('flow');

psthN = mean(N(:,:,clu),2);
psthTaskFit = mean(task_fit(:,:,clu),2);
for i = 1:numel(trix)
    trialdatN(:,i) = mySmooth(N(:,trix(i),clu),sm)./5;
    trialdatTaskFit(:,i) = mySmooth(task_fit(:,trix(i),clu),sm);
end

ax = nexttile; hold on;
patchline(obj.time,psthN,'EdgeColor',[0 0 0],'LineWidth',2,'EdgeAlpha',0.3);
patchline(obj.time,psthTaskFit,'EdgeColor',[41, 128, 64]./255,'LineWidth',2,'EdgeAlpha',1)
xline(obj.time(goCueIx),'k--','LineWidth',gclw)
xlim([obj.time(10) obj.time(end)])
title(['PSTH ' num2str(clu)])
ax.FontSize = 13;

for i = 1:numel(trix)
    ax = nexttile; hold on;
    patchline(obj.time,trialdatN(:,i),'EdgeColor',[0 0 0],'LineWidth',2,'EdgeAlpha',0.3);
    patchline(obj.time,trialdatTaskFit(:,i),'EdgeColor',[41, 128, 64]./255,'LineWidth',1,'EdgeAlpha',1)
    xline(obj.time(goCueIx),'k--','LineWidth',gclw)
    xlim([obj.time(10) obj.time(end)])
    title(['Trial ' num2str(trix(i))])
    ax.FontSize = 13;
end



title(t, 'task fit reconstructions')
xlabel(t, 'Time (s) from go cue')
ylabel(t, 'Centered Firing Rate (Hz)')

if sav
    pth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\musall2019\figs';
    fn = ['clu' num2str(clu) '_taskfit_recon'];
    mysavefig(f1,pth,fn)
end

% reconstruct PSTHs and few single trials, ui fit

% N is original data in (time*trials,neurons)
% fullFit is reocnstructed data of same size

ui_fit = reshape(ui_fit,sz.N(1),sz.N(2),sz.N(3));

%


f2 = figure;
f2.Position = [867    80   379   897];
t = tiledlayout('flow');

psthUIFit = mean(ui_fit(:,:,clu),2);
for i = 1:numel(trix)
    trialdatUIFit(:,i) = mySmooth(ui_fit(:,trix(i),clu),sm);
end

ax = nexttile; hold on;
patchline(obj.time,psthN,'EdgeColor',[0 0 0],'LineWidth',2,'EdgeAlpha',0.3);
patchline(obj.time,psthUIFit,'EdgeColor',[214, 47, 125]./255,'LineWidth',2,'EdgeAlpha',1)
xline(obj.time(goCueIx),'k--','LineWidth',gclw)
xlim([obj.time(10) obj.time(end)])
title(['PSTH ' num2str(clu)])
ax.FontSize = 13;

for i = 1:numel(trix)
    ax = nexttile; hold on;
    patchline(obj.time,trialdatN(:,i),'EdgeColor',[0 0 0],'LineWidth',2,'EdgeAlpha',0.3);
    patchline(obj.time,trialdatUIFit(:,i),'EdgeColor',[214, 47, 125]./255,'LineWidth',1,'EdgeAlpha',1)
    xline(obj.time(goCueIx),'k--','LineWidth',gclw)
    xlim([obj.time(10) obj.time(end)])
    title(['Trial ' num2str(trix(i))])
    ax.FontSize = 13;
end



title(t, 'Unin. Move fit reconstructions')
xlabel(t, 'Time (s) from go cue')
ylabel(t, 'Centered Firing Rate (Hz)')

if sav
    pth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\musall2019\figs';
    fn = ['clu' num2str(clu) '_uifit_recon'];
    mysavefig(f2,pth,fn)
end

% reconstruct PSTHs and few single trials, i fit

% N is original data in (time*trials,neurons)
% fullFit is reocnstructed data of same size

i_fit = reshape(i_fit,sz.N(1),sz.N(2),sz.N(3));

%


f3 = figure;
f3.Position = [1250          74         379         897];
t = tiledlayout('flow');

psthIFit = mean(i_fit(:,:,clu),2);
for i = 1:numel(trix)
    trialdatIFit(:,i) = mySmooth(i_fit(:,trix(i),clu),sm);
end

ax = nexttile; hold on;
patchline(obj.time,psthN,'EdgeColor',[0 0 0],'LineWidth',2,'EdgeAlpha',0.3);
patchline(obj.time,psthIFit,'EdgeColor',[50, 62, 219]./255,'LineWidth',2,'EdgeAlpha',1)
xline(obj.time(goCueIx),'k--','LineWidth',gclw)
xlim([obj.time(10) obj.time(end)])
title(['PSTH ' num2str(clu)])
ax.FontSize = 13;

for i = 1:numel(trix)
    ax = nexttile; hold on;
    patchline(obj.time,trialdatN(:,i),'EdgeColor',[0 0 0],'LineWidth',2,'EdgeAlpha',0.3);
    patchline(obj.time,trialdatIFit(:,i),'EdgeColor',[50, 62, 219]./255,'LineWidth',1,'EdgeAlpha',1)
    xline(obj.time(goCueIx),'k--','LineWidth',gclw)
    xlim([obj.time(10) obj.time(end)])
    title(['Trial ' num2str(trix(i))])
    ax.FontSize = 13;
end



title(t, 'In. Move fit reconstructions')
xlabel(t, 'Time (s) from go cue')
ylabel(t, 'Centered Firing Rate (Hz)')

if sav
    pth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements\musall2019\figs';
    fn = ['clu' num2str(clu) '_ifit_recon'];
    mysavefig(f3,pth,fn)
end


%% coding directions

clear temp

% task
dat = task_fit; % (time,trials,clu)
temp{1} = squeeze(mean(dat(:,trials.rightix,:),2));
temp{2} = squeeze(mean(dat(:,trials.leftix,:),2));
input_data = cat(3,temp{1},temp{2});
getCDs(obj,input_data,params,'Task')

% instructed
dat = i_fit; % (time,trials,clu)
temp{1} = squeeze(mean(dat(:,trials.rightix,:),2));
temp{2} = squeeze(mean(dat(:,trials.leftix,:),2));
input_data = cat(3,temp{1},temp{2});
getCDs(obj,input_data,params,'Instructed Movements')

% uninstructed
dat = ui_fit; % (time,trials,clu)
temp{1} = squeeze(mean(dat(:,trials.rightix,:),2));
temp{2} = squeeze(mean(dat(:,trials.leftix,:),2));
input_data = cat(3,temp{1},temp{2});
getCDs(obj,input_data,params,'Uninstructed Movements')








