% Predicting context from populations of single cells

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
addpath(genpath(fullfile(utilspth,'fig3')));
figpth = [basepth '\Uninstructed-Movements\Fig 3'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context switching')));
figpth = [basepth '\Uninstructed-Movements\Fig 6'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context_funcs')));
figpth = [basepth '\Uninstructed-Movements\Fig 5'];
addpath(genpath(fullfile(figpth,'funcs')));

load([basepth '\Uninstructed-Movements\ContextColormap.mat']);
%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials

params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};                   % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};                   % error left, no stim, aw off
params.condition(end+1) = {'~early&~stim.enable&~autowater'};                          %  no stim, aw off

params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};             % right hits, no stim, aw on
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};             % left hits, no stim, aw on
params.condition(end+1) = {'R&miss&~stim.enable&autowater'};                   % error right, no stim, aw on
params.condition(end+1) = {'L&miss&~stim.enable&autowater'};                   % error left, no stim, aw on
params.condition(end+1) = {'~early&~stim.enable&autowater'};                          %  no stim, aw on

params.tmin = -3;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

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
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);

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
%% context decoding from population of single neurons

clearvars -except datapth kin me meta obj params acc acc_shuf

% params
rez.nFolds = 4; % number of iterations (bootstrap)

rez.binSize = 50; % ms
rez.dt = floor(rez.binSize / (params(1).dt*1000)); % samples
rez.tm = obj(1).time(1:rez.dt:numel(obj(1).time));
rez.numT = numel(rez.tm);

rez.train = 1; % fraction of trials to use for training (1-train for testing)

rez.nShuffles = 20;

% match number of right and left 2afc, and right and left aw
cond2use = 2:11;
afccond = [1 2];
awcond = [6 7];

for sessix = 1:numel(obj)
    disp(['Decoding session ' num2str(sessix) ' / ' num2str(numel(obj))])

    % trials
    trials_cond = params(sessix).trialid(cond2use);

    minAFCTrials = cellfun(@(x) numel(x),trials_cond(afccond), 'UniformOutput',false);
    nafc = min(cell2mat(minAFCTrials));

    minAWTrials = cellfun(@(x) numel(x),trials_cond(awcond), 'UniformOutput',false);
    naw = min(cell2mat(minAWTrials));

    nTrials = min(nafc,naw);

    trials_afc = cellfun(@(x) randsample(x,nTrials), trials_cond(afccond), 'UniformOutput', false);
    trialsAFC = cell2mat(trials_afc);
    trialsAFC = trialsAFC(:);

    trials_aw = cellfun(@(x) randsample(x,nTrials), trials_cond(awcond), 'UniformOutput', false);
    trialsAW = cell2mat(trials_aw);
    trialsAW = trialsAW(:);

    trials.all = [trialsAFC ; trialsAW];

    % labels (1 for afc, -1 for aw)
    Y = [ones(nTrials,1) ; ones(nTrials,1) ; -ones(nTrials,1) ; -ones(nTrials,1)]; % right afc, left afc, left aw, right aw
    %         Y = [ones(nafc,1) ; ones(nafc,1) ; -ones(naw,1) ; -ones(naw,1)]; % right hits, left hits, left miss, right miss

    % input
    % get single-trial PSTHs
    temp = obj(sessix).trialdat(:,:,trials.all);                       % (time x cells x trials)
    X = permute(temp,[1 3 2]);                                         % Re-order to (time x trials x cells)

    % train/test split
    [trials.train,trials.trainidx] = datasample(trials.all,round(numel(trials.all)*rez.train),'Replace',false);
    trials.testidx = find(~ismember(trials.all,trials.train));
    trials.test = trials.all(trials.testidx);

    in.train.y = Y(trials.trainidx);
    in.test.y  = Y(trials.testidx);
    in.train.X = X(:,trials.trainidx,:);
    in.test.X  = X(:,trials.testidx,:);

    % decoding
    acc(:,sessix) = NeuronContext_Decoder(in,rez,trials);


    % shuffle labels for a 'null' distribution
    Y = randsample(Y,numel(Y));

    % train/test split
    in.train.y = Y(trials.trainidx);
    in.test.y  = Y(trials.testidx);

    for ishuf = 1:rez.nShuffles
        acc_shuf(:,sessix,ishuf) = NeuronContext_Decoder(in,rez,trials);
    end
end

acc_shuf_ = reshape(acc_shuf,size(acc_shuf,1),size(acc_shuf,2)*size(acc_shuf,3));
%% plot

cols = {'k',[0.6,0.6,0.6]};

alph = 0.15;

sample = mode(obj(1).bp.ev.sample) - mode(obj(1).bp.ev.goCue);
delay = mode(obj(1).bp.ev.delay) - mode(obj(1).bp.ev.goCue);
trialStart = mode(obj(1).bp.ev.bitStart) - mode(obj(1).bp.ev.goCue);

figure;
ax = gca;
hold on;
smooth = 21;
toplot = mySmooth(mean(acc,2),smooth);
err = 1.96*mySmooth((std(acc,[],2)./sqrt(numel(obj))),smooth);
shadedErrorBar(rez.tm(1:end-1),toplot,err,{'Color',cols{1},'LineWidth',2},alph,ax)

toplot = mySmooth(mean(acc_shuf_,2),smooth);
err = 1.96*mySmooth((std(acc_shuf_,[],2)./sqrt(numel(obj))),smooth);
shadedErrorBar(rez.tm(1:end-1),toplot,err,{'Color',cols{2},'LineWidth',2},alph,ax)
% shadedErrorBar(rez.tm(1:end-1),mean(acc,2),getCI(acc),{'Color',cols{1},'LineWidth',2},alph,ax)
% shadedErrorBar(rez.tm(1:end-1),mean(acc_shuf_,2),getCI(acc_shuf_),{'Color',cols{2},'LineWidth',2},alph,ax)
xline(0,'k--','LineWidth',1)
xline(sample,'k--','LineWidth',1)
xline(delay,'k--','LineWidth',1)
xlim([trialStart, params(1).tmax-0.2])
ylim([ax.YLim(1) 1])

xlabel('Time from go cue (s)')
ylabel([num2str(rez.nFolds) '-Fold CV Accuracy'])
title('Context Decoding from Population Activity')

h = zeros(2, 1);
for i = 1:numel(h)
    h(i) = plot(NaN,NaN,'-','Color',cols{i},'LineWidth',2);
end
legString = {'ctrl','shuf'};

leg = legend(h, legString);
leg.EdgeColor = 'none';
leg.Location = 'best';

ax.FontSize = 15;

