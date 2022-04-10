clear,clc,close all

addpath(genpath(pwd))


%% SET RUN PARAMS
params.alignEvent          = 'jawOnset'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 0.5; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};          % right hits, no stim, aw on
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};          % left hits, no stim, aw on
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};        % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};        % error left, no stim, aw off
params.condition(end+1) = {'~hit&~miss&~stim.enable&~autowater&~early'};    % ignore
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};           % hit 2afc
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};            % hit aw


% set conditions used for finding activity modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %rhit, aw off 
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %lhit, aw off 
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %rmiss, aw off 
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %lmiss, aw off 
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};    % hit, aw off 

% specify probe number and areas to load and process data
params.probe(1) = 1;
params.probeArea{1} = 'ALM';

% params.probe(end+1) = 2; % for multiprobe stuff, not ready yet (TODO)
% params.probeArea{end+1} = 'ALM';

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

params.vid_dt = 1/400;

%% SET METADATA
% experiment meta data
meta.datapth = '/Users/Munib/Documents/Economo-Lab/code/data';
% meta.datapth = '/Volumes/MUNIB_SSD/Experiments';
meta.anm = 'EKH3'; % 'JEB7'  'EKH3'
meta.date = '2021-08-11'; % '2021-04-29'   '2021-08-11
meta.datafn = findDataFn(meta);

%% LOAD AND PROCESS DATA
dat = load(fullfile(meta.datapth,meta.datafn));
obj = dat.obj;

for prbnum = 1:numel(params.probe)
    disp('______________________________________________________')
    disp(['Processing data for probe ' num2str(prbnum)])
    disp(' ')
    [params,obj] = processData(obj,params,prbnum);
end

% if only one probe, clean up so all previous code works
if numel(params.probe)==1
    [obj,params] = oneProbeTrim(obj,params);
end

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')

selectivity = obj.psth(:,:,1) - obj.psth(:,:,2);
figure; imagesc(selectivity')
colorbar
caxis([-40,40])

%% label move or non-move
% [obj.movix,obj.movtime] = getMoveIdx(obj,params);

%% ACTIVITY MODES
rez.time = obj.time;
rez.psth = obj.psth;
rez.condition = params.condition;
rez.alignEvent = params.alignEvent;

%% jaw mode
fr = cat(3, obj.trialdat(:,:,params.trialid{1}), obj.trialdat(:,:,params.trialid{2}));
fr(isnan(fr)) = 0;


traj = obj.traj{1};  % Get the video data

conds = [1 2];
trials = {params.trialid{conds}};
trials = cell2mat(trials');
taxis = params.tmin:params.dt:params.tmax;

jaw = nan(numel(obj.time),numel(trials));
for i = 1:numel(trials)                        % For every trial in the condition
    trix = trials(i);
    if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
        continue;
    end
    
    if ~isnan(traj(trix).frameTimes)                           % If the video data from this trial is good...
        ts = mySmooth(traj(trix).ts(:, 2, 2), 21);                                               % Side-view, up and down position of the jaw, smoothed
        tsinterp = interp1(traj(trix).frameTimes-0.5-(obj.bp.ev.jawOnset(trix)), ts, taxis);          % Linear interpolation of jaw position to keep number of time points consistent across trials
        basederiv = median(myDiff(tsinterp,400),'omitnan');                                         % Find the median jaw velocity (aka baseline)
    end
    %Find the difference between the jaw velocity and the
    %baseline jaw velocity
    jaw(:, i) = abs(diff(tsinterp)-basederiv);      % Values > 0 = jaw is moving
end



rez.jaw_mode = zeros(size(fr, 2), 1);
for i = 1:size(fr, 2)
    f = mySmooth(squeeze(fr(:, i, :)), 51);
    j = jaw;
    j(isnan(j)) = 0;
    tmp = corrcoef(f(:), j(:));
    rez.jaw_mode(i) = tmp(1,2);
end

%% choice mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
cond{3} = params.modecondition{3};
cond{4} = params.modecondition{4};
epoch = 'delay';
rez.choice_mode = choiceMode(obj,params,cond,epoch,rez.alignEvent);
clear cond

%% remainder mode
orthModes = [rez.jaw_mode rez.choice_mode];
modesToKeep = eye(size(obj.psth,2)) - (orthModes*orthModes');

residualpsth = nan(size(obj.psth));
for i = 1:size(obj.psth,3)
    residualpsth(:,:,i) = obj.psth(:,:,i) * modesToKeep;
end

X = [residualpsth(:,:,1) ; residualpsth(:,:,2)]; % left and right 2afc

% SVD
nComp = 6;
V = myPCA(X - mean(X));
for i = 1:nComp
    rez.(['remainder' num2str(i) '_mode']) = V(:,i); % decreasing order of var explained
end


%% orthogonalize

[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
modes = zeros(numel(params.cluid),numel(fns));
for i = 1:numel(fns)
    modes(:,i) = rez.(fns{i});
end

orthModes = gschmidt(modes);

for i = 1:numel(fns)
    rez.(fns{i}) = orthModes(:,i);
end


%% PLOT MODES

% MODES VIZ

% plot correct trials and AW trials
plt.title = '';
plt.legend = {'Right 2AFC','Left 2AFC'};
plt.conditions = [1,2];
plt.lw = [2.7 2.7];
plt.smooth = 31;
plt.colors = {[0 0 1],[1 0 0]};
plt.save = 0;
plotAllModes(rez, obj.bp.ev, params.alignEvent, plt) 

% ORTHOGONALITY VIZ
% dotProductModes(rez,modes,'NOT ORTHOGONALIZED')
% dotProductModes(rez,orthModes,'ORTHOGONALIZED')


%% get trial-averaged latents and single-trial latents

% trial-avg latents
% latents_avg is a (time,nCond*nModes) input matrix

[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
conds = [1 2];
latents_avg = zeros(size(rez.psth,1),numel(fns)*numel(conds));
ct = 1;
for i = 1:numel(fns)
    for j = 1:numel(conds)
        cond = conds(j);
        latents_avg(:,ct) = rez.psth(:,:,cond)*rez.(fns{i});
        ct = ct + 1;
    end
end

% single trial latents
% latents is a (time,nModes,trials) input matrix
sm = 101;
trials = {params.trialid{conds}};
trials = cell2mat(trials');
latents = zeros(size(rez.psth,1),numel(fns),numel(trials));
for i = 1:numel(fns)
    for j = 1:numel(trials)
        latents(:,i,j) = mySmooth(obj.trialdat(:,:,trials(j)) * rez.(fns{i}),sm);
    end
end
% close all; figure; imagesc(squeeze(latents(:,3,:))'); colorbar; 

%%
% 
% % output matrix (to be predicted) is a bunch of features (motion energy,dlc
% % feature position, velocity, etc.)
% featNames = getFeatNames(obj);
% 
% feat2use = {[2],
%              []}; % which features to use for each view. If empty array, won't use that view
% 
% % feature trajectories
% traj = getFeatTraj(obj,featNames,feat2use);
% 
% % align feature trajectories to move onset (estimated with jaw position)
% % trials that don't have t_preEv or t_postEv time points to use are
% % discarded. See trialid field for the actual trial number. 
% traj = alignFeatTraj(traj,obj.bp.ev.jawOnset,abs(params.tmin),params.tmax,params.vid_dt);
% 
% % trim to trials from previous section only
% discard = [];
% for i = 1:numel(traj)
%     if ~ismember(traj(i).trialid,trials)
%         discard = [discard ; i];
%     end
% end
% traj(discard) = [];
% 
% % pull out jaw pos into a matrix
% jawPos = zeros(numel(traj(1).time),2,numel(traj)); % (time,nFeat*nCoord,trials)
% for trix = 1:numel(traj)
%     jawPos(:,:,trix) = [traj(trix).jaw_cam0(:,1:2)];
% end
% 
% % fillnans
% jawPos = fillmissing(jawPos,'nearest');
% 
% % downsample jaw pos to match clu (viddt = 400, cludt = 200)
% jawPos_downsamp = zeros(ceil(size(jawPos,1)/2),size(jawPos,2),size(jawPos,3));
% for trix = 1:numel(traj)
%     jawPos_downsamp(:,:,trix) = downsample(jawPos(:,:,trix),2);
% end
% jawPos = jawPos_downsamp(1:end-1,:,:);
% 
% % jawVel
% dt = params.vid_dt;
% sm = 25;
% jawVel = zeros(size(jawPos));
% for trix = 1:size(jawPos,3)
%     for cix = 1:size(jawPos,2)
%         temp = myDiff(jawPos(:,cix,trix),dt);
%         jawVel(:,cix,trix) = sgolayfilt(temp,3,sm);
%     end
% end
% 
% % jawAcc
% dt = params.vid_dt;
% sm = 25;
% jawAcc = zeros(size(jawPos));
% for trix = 1:size(jawPos,3)
%     for cix = 1:size(jawPos,2)
%         temp = myDiff(jawVel(:,cix,trix),dt);
%         jawAcc(:,cix,trix) = sgolayfilt(temp,3,sm);
%     end
% end
% 
% % close all; figure; imagesc(squeeze(jawAcc(:,2,:))'); colorbar; %caxis([-1000 1000])
% 
% % concatenate time and trials to make each 3d matrix a 2d matrix of size
% % (time*trials,numFeats)
% latents = permute(latents,[1 3 2]);
% latents = reshape(latents,size(latents,1)*size(latents,2),size(latents,3));
% 
% jawPos = permute(jawPos,[1 3 2]);
% jawPos = reshape(jawPos,size(jawPos,1)*size(jawPos,2),2);
% 
% jawVel = permute(jawVel,[1 3 2]);
% jawVel = reshape(jawVel,size(jawVel,1)*size(jawVel,2),2);
% 
% jawAcc = permute(jawAcc,[1 3 2]);
% jawAcc = reshape(jawAcc,size(jawAcc,1)*size(jawAcc,2),2);
% 
% 
% save('jeb7_decoding','latents_avg','latents','jawPos','jawVel','jawAcc')
% 
% 
% %% 
% 
% % feature velocities
% sm = 25; % smoothing window for sgolay filter
% vel = getFeatVelocity(traj,params.vid_dt,sm);
% 
% 
% % create a feature matrix
% featmat = zeros(numel(traj(1).time),numel(fieldnames(vel))*2,numel(vel)); % (time,nFeat*nCoord,trials)
% for trix = 1:numel(vel)
%     [featmat(:,:,trix),labels] = featStruct2Mat(vel(trix));
% end
% 
% 
% 
% 
% %% morlet wavelet transform
% % 
% % %%%%%%%% Wavelet Parameters %%%%%%%%
% % %number of wavelet frequencies to use
% % numPeriods = 25;
% % %dimensionless Morlet wavelet parameter
% % omega0 = 5;
% % %minimum frequency for wavelet transform (Hz)
% % minF = 0.5;
% % %maximum frequency for wavelet transform (Hz)
% % maxF = 60;
% % 
% % minT = 1 ./ maxF;
% % maxT = 1 ./ minF;
% % Ts = minT.*2.^((0:numPeriods-1).*log(maxT/minT)/(log(2)*(numPeriods-1)));
% % f = fliplr(1./Ts);
% % 
% % amplitudes = struct();
% % fnames = fieldnames(vel);
% % for trix = 1:numel(traj)
% %     for featix = 1:numel(fnames,2)
% % 
% %         amplitudes(1).(v
% %         amplitudes{featix} = fastWavelet_morlet_convolution_parallel(featmat(:,featix,trix),f,omega0,params.vid_dt);
% %     end
% % end
% %%
% 
% % downsample feature matrix (viddt = 400, cludt = 200)
% jawPos_downsamp = zeros(ceil(size(featmat,1)/2),size(featmat,2),size(featmat,3));
% for trix = 1:numel(vel)
%     jawPos_downsamp(:,:,trix) = downsample(featmat(:,:,trix),2);
% end
% featmat = jawPos_downsamp(1:end-1,:,:);
% 
% %%
% % % decoding examples
% % ridgeParam = 70000; % need a different value for each feature!!!!!
% % trix = 1;
% % figure; hold on
% % for i = 1:size(featmat,2)
% %     W = ridge(testfeat,latents,ridgeParam);
% %     testfeat = envelope(fillmissing(featmat(:,i,trix),'nearest'),200,'rms'); % upper envelope via hilbert transform
% %     testfeat_ = latents*W;
% %     a = 0;
% %     b = 1;
% %     col = (b-a).*rand(3,1) + a;
% %     plot(testfeat,'Color',col,'LineWidth',2);
% %     plot(testfeat_,'Color',col,'LineWidth',2)
% %     
% % end
% % 
% % figure; 
% % k = 70000;
% % clear error
% % for i = 1:numel(k)
% %     W = ridge(testfeat,latents,k(i));
% %     testfeat_ = latents*W;
% %     error(i) = immse(testfeat,testfeat_);
% %     
% % end
% % figure; plot(k,error);
% % 
% % 
% % 
% % figure;
% % for i = 100:100:1000
% %     testfeat = envelope(featmat(:,1,1),i,'rms');
% %     plot(testfeat)
% %     title(num2str(i))
% %     pause
% %     clf
% % end



