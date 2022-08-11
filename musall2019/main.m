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
params.reload = false; % redo video svd

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
params.trials2use = params.trialid{3}; % right and left, hit,miss, ~stim,~autowater

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')


%%
clearvars -except meta params obj datapth

%% load raw video data if needed


vidpth = fullfile(datapth,'Video',meta.anm,meta.date,'Cam0');
vids = dir(vidpth);
vids = {vids.name}';

vids = patternMatchCellArray(vids,{meta.anm,meta.date},'all');
vids = natsortfiles(vids);

vids = vids(params.trials2use); % trim trials

vidtm = params.tmin:params.viddt:params.tmax;
nFrames = numel(vidtm);

fpth = fullfile(datapth,'Video',meta.anm,meta.date);
if ~exist(fullfile(fpth ,'vidR.mat'), 'file') || ~exist(fullfile(fpth ,'absVidR.mat'), 'file') || params.reload

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
    vidSize = size(grayframes);
    grayframes = reshape(grayframes,vidSize(1),vidSize(2),vidSize(3)*vidSize(4)); % (x,y,frames*trials)
end



%% check for video SVD and create if missing

nrDims = 200; % svd dims to keep
if ~exist(fullfile(vidpth, 'vidR.mat'), 'file') || params.reload
    cVideo = reshape(grayframes, [], size(grayframes,3)); %merge all pixels

    [~,s,Vr] = fsvd(single(cVideo),nrDims,1,0,0);
    vidR = s*Vr';
    vidR = bsxfun(@minus, vidR, mean(vidR, 2));
    vidR = bsxfun(@rdivide, vidR, std(vidR, [], 2))'; %this is the video regressor in the model
    vidR = reshape(vidR,vidSize(3),vidSize(4),nrDims); % (time,trials,videodims)
    save(fullfile(vidpth, 'vidR.mat'),'vidR');
else
    load(fullfile(vidpth ,'vidR.mat'),'vidR');
end


%% compute motion energy and repeat svd if needed
if ~exist(fullfile(vidpth ,'absVidR.mat'), 'file') || params.reload
    cVideo = grayframes; %merge all pixels
    cVideo = cat(3,cVideo(:,:,1),cVideo); %duplicate first frame
    cVideo = abs(diff(cVideo,[],3)); %compute absolute of temporal derivative (motion energy)
    cVideo = reshape(cVideo, [], size(cVideo,3)); %merge all pixels

    [~, s, Vr] = fsvd(single(cVideo), nrDims, 1, 0, 0);
    absVidR = s * Vr';
    absVidR = bsxfun(@minus, absVidR, mean(absVidR, 2));
    absVidR = bsxfun(@rdivide, absVidR, std(absVidR, [], 2))'; %this is the motion energy regressor in the model
    absVidR = reshape(absVidR,vidSize(3),vidSize(4),nrDims); % (time,trials,videodims)
    save(fullfile(vidpth, 'absVidR.mat'),'absVidR');
else
    load(fullfile(vidpth ,'absVidR.mat'),'absVidR');
end
clear cVideo


%% synchronize neural and video data

clutm = obj.time;

vidtm = params.tmin:params.viddt:params.tmax;

vidR = syncTimes(vidtm,vidR,clutm);

absVidR = syncTimes(vidtm,absVidR,clutm);


%% ephys data

N = permute(obj.trialdat,[1 3 2]); % single trial binned neural data (time,trials,clusters)

% trim trials
N = N(:,params.trials2use,:);

sz.N = size(N);
sz.VidR = size(vidR);
sz.AbsVidR = size(absVidR);
sz.leg = {'time','trials','clu/dims'};


%% reshape video and ephys data

vidR = reshape(vidR,size(vidR,1)*size(vidR,2),size(vidR,3)); % (time*trials,viddims)

absVidR = reshape(absVidR,size(absVidR,1)*size(absVidR,2),size(absVidR,3)); % (time*trials,viddims)

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

% uninstructed licks
uiLicks = lick;
[~,goCueIx] = min(abs(obj.time - 0));
uiLicks(goCueIx:end,:) = 0;
moveEvents(:,1) = uiLicks(:);

% instructed licks
iLicks = lick;
iLicks(1:(goCueIx-1),:) = 0;
moveEvents(:,2) = iLicks(:);


%% make movement design matrix

moveLabels = {'preGoCueLicks' 'postGoCueLicks'}; % some movement variables
moveEventType = [3 3]; %different type of events.

[moveR, moveIdx] = makeDesignMatrix(moveEvents, moveEventType, params); %make design matrix for task variables


% figure; plot(moveR(:,moveIdx==1),'r')
% hold on
% plot(moveR(:,moveIdx==2),'g')

%% video regressors

% make video traces before go cue uninstructed and after go cue instructed



%% make full design matrix

fullR = [taskR, moveR, absVidR, vidR]; %make new, single design matrix

fullR = bsxfun(@minus, fullR, mean(fullR, 1));

moveLabels = [moveLabels, {'ME'} , {'video'}];

% labels for different regressor sets. It is REALLY important this agrees with the order of regressors in fullR.
regLabels = cat(2,taskLabels,moveLabels);

%index to reconstruct different response kernels
regIdx = [
            ones(sum(taskIdx==1),1) * find(ismember(regLabels,'choice')); ...
            ones(sum(taskIdx==2),1) * find(ismember(regLabels,'prevchoice')); ...
            ones(sum(taskIdx==3),1) * find(ismember(regLabels,'success')); ...
            ones(sum(taskIdx==4),1) * find(ismember(regLabels,'prevsuccess')); ...
            ones(sum(taskIdx==5),1) * find(ismember(regLabels,'reward')); ...
            ones(sum(moveIdx==1),1) * find(ismember(regLabels,'preGoCueLicks')); ...
            ones(sum(moveIdx==2),1) * find(ismember(regLabels,'postGoCueLicks')); ...
            ones(nrDims,1) * find(ismember(regLabels,'ME')); ...
            ones(nrDims,1) * find(ismember(regLabels,'video')) ...
         ];


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
stimIdx = regIdx == find(ismember(regLabels,'stim'));

uninstructedLabels = {'preGoCueLicks'};
lbIdx = find(ismember(regLabels,uninstructedLabels));
uiIdx = ismember(regIdx ,lbIdx) ; % logical of uninstructed vars in cols of fullR

not_ui_fit = fullR(:, ~uiIdx) * dimBeta(~uiIdx, :); %fit data
ui_fit = fullR(:, uiIdx) * dimBeta(uiIdx, :); %fit data


disp('Full model fit completed');


%% variance explained (https://economictheoryblog.com/2014/11/05/proof/)
R = corrcoef(N,fullFit);
ve_full = R(1,2)^2;

R = corrcoef(N,not_ui_fit);
ve_not_ui = R(1,2)^2;

R = corrcoef(N,ui_fit);
ve_ui = R(1,2)^2;

figure;
X = categorical({'Full','Not UI','UI'});
bar(X, [ve_full,ve_not_ui,ve_ui])

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




