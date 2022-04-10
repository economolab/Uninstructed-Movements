clear,clc,close all

addpath(genpath(pwd))

% load lfads input and output
datapth = '/Users/Munib/Documents/Economo-Lab/code/lfads/input_data/';
datafn = 'JEB7_2021-04-29_run2.mat';
temp = load(fullfile(datapth,datafn));
datin = temp.lfads; clear temp;

% load lfads output
lfadsdatapth = '/Users/Munib/Documents/Economo-Lab/code/lfads/output';
lfadsfname = 'model_runs_h5_train_posterior_sample_and_average_JEB7-2021-04-29_run2';
datout = HDF5ToStruct(lfadsdatapth,lfadsfname);

% conditions
rhit = find(datin.obj.bp.R & datin.obj.bp.hit & ~datin.obj.bp.early & ~datin.obj.bp.autowater);
lhit = find(datin.obj.bp.L & datin.obj.bp.hit & ~datin.obj.bp.early & ~datin.obj.bp.autowater);
rmiss = find(datin.obj.bp.R & datin.obj.bp.miss & ~datin.obj.bp.early & ~datin.obj.bp.autowater);
lmiss = find(datin.obj.bp.L & datin.obj.bp.miss & ~datin.obj.bp.early & ~datin.obj.bp.autowater);
if contains(lfadsfname,'train')
    trials = datin.train_trials;
    
elseif contains(lfadsfname,'valid')
    trials = datin.valid_trials;
end
rhit = ismember(trials,rhit);
lhit = ismember(trials,lhit);
rmiss = ismember(trials,rmiss);
lmiss = ismember(trials,lmiss);

%% ACTIVITY MODES
rez.time = datin.obj.time;
rez.psth = datin.obj.psth;
rez.condition = datin.preprocess_params.condition;
rez.alignEvent = datin.preprocess_params.alignEvent;

%% jaw mode
fr = datout.output_dist_params; % fr from lfads
% fr = datout.factors; % single trial projections onto factors from lfads
fr(isnan(fr)) = 0;
fr = permute(fr,[2,1,3]);

traj = datin.obj.traj{1};  % Get the video data

trials = datin.train_trials;
taxis = datin.obj.time;

jaw = nan(numel(taxis),numel(trials));
for i = 1:numel(trials)                        % For every trial in the condition
    trix = trials(i);
    if isnan(traj(trix).NdroppedFrames )                       % If the video data from this trial isn't good, skip it
        continue;
    end
    
    if ~isnan(traj(trix).frameTimes)                           % If the video data from this trial is good...
        ts = mySmooth(traj(trix).ts(:, 2, 2), 21);                                               % Side-view, up and down position of the jaw, smoothed
        tsinterp = interp1(traj(trix).frameTimes-0.5-(datin.obj.bp.ev.jawOnset(trix)), ts, taxis);          % Linear interpolation of jaw position to keep number of time points consistent across trials
        tsinterp = fillmissing(tsinterp,'nearest');
        basederiv = median(myDiff(tsinterp,400),'omitnan');                                         % Find the median jaw velocity (aka baseline)
    end
    %Find the difference between the jaw velocity and the
    %baseline jaw velocity
    jaw(:, i) = abs(myDiff(tsinterp,400)-basederiv);      % Values > 0 = jaw is moving
end



rez.jaw_mode = zeros(size(fr, 2), 1);
for i = 1:size(fr, 2) % for each neuron
    f = mySmooth(squeeze(fr(:, i, :)), 51);
    j = jaw;
    j(isnan(j)) = 0;
    tmp = corrcoef(f(:), j(:));
    rez.jaw_mode(i) = tmp(1,2);
end

jaw_latent = zeros(size(fr,1),size(fr,3));
for i = 1:size(fr,3)
    jaw_latent(:,i) = squeeze(fr(:,:,i)) * rez.jaw_mode;
end

trialOrder = [find(rhit);find(lhit);find(rmiss);find(lmiss)];
jaw = jaw(:,trialOrder);
jaw_latent = jaw_latent(:,trialOrder);

figure(1); clf(figure(1));
ax(1) = subplot(121);
imagesc(datin.obj.time,1:size(fr,2),jaw')
colorbar(ax(1))
caxis(ax(1),[0 0.03]);
ax(2) = subplot(122);
imagesc(datin.obj.time,1:size(fr,2),jaw_latent')
colorbar(ax(2))
% caxis(ax(2),[0 0.06]);
linkaxes(ax,'xy');
cc = corrcoef(jaw(:),jaw_latent(:));
sgtitle(num2str(cc(1,2)))


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

