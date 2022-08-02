clear, clc, close all

addpath(genpath(pwd))

%% TODO
% - created new objs from pipeline and nidq data - only use those
% - need to handle frameTimes appropriately now (subtract 0.5 secs)
% - for choice decoding, use ~50ms bins, rather than single time points

%% default params (same for each session)

dfparams = [];

% dfparams.alignEv = 'goCue';
% dfparams.times = [-2.5 1.5]; % relative to goCue
dfparams.alignEv = 'delay';
dfparams.times = [-1.0 2.5]; % relative to delay

dfparams.dt_vid = 0.0025;
dfparams.time = dfparams.times(1):dfparams.dt_vid:dfparams.times(2);

dfparams.cond(1) = {'(hit|miss|no)&~stim.enable&~autowater&~autolearn'}; % right and left trials, no stim, no autowater
dfparams.cond(end+1) = {'(hit|miss|no)&stim.enable&~autowater&~autolearn'};  % right and left trials, stim, no autowater
dfparams.cond(end+1) = {'R&~stim.enable&~autowater&~autolearn'}; % right trials, no stim, no autowater
dfparams.cond(end+1) = {'R&stim.enable&~autowater&~autolearn'};  % right trials, stim, no autowater
dfparams.cond(end+1) = {'L&~stim.enable&~autowater&~autolearn'}; % left trials, no stim, no autowater
dfparams.cond(end+1) = {'L&stim.enable&~autowater&~autolearn'};  % left trials, stim, no autowater

dfparams.plt.color{1}     = [10, 10, 10];
dfparams.plt.color{end+1} = [120, 117, 117];
dfparams.plt.color{end+1} = [31, 23, 252];
dfparams.plt.color{end+1} = [22, 172, 247];
dfparams.plt.color{end+1} = [252, 23, 23];
dfparams.plt.color{end+1} = [252, 23, 130];
dfparams.plt.color = cellfun(@(x) x./255, dfparams.plt.color, 'UniformOutput', false);
dfparams.plt.ms = {'.','.','x','x','o','o'};

%% load data objects

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];
meta = loadMAH13_MCStim(meta,datapth);


obj = loadObjs(meta);

%% find trials for each condition

for i = 1:numel(meta)
    params(i).trialid = findTrials(obj(i), dfparams.cond);
end

%% behavioral performance

rez = getPerformance(meta,obj,params);

% plots
plotPerformanceAllMice(meta,obj,rez,dfparams,params)
% plotPerformanceEachMouse(meta,obj,rez,dfparams,params) % TODO


%% number of early licks per trial

% TODO: classify early licks by kinematics, not bpod data

cond2use = 3:6;
early_per_trial = nan(numel(obj),numel(cond2use)); % (sessions,conds)
for sessix = 1:numel(obj)
    for condix = 1:numel(cond2use)
        trials2use = params(sessix).trialid{cond2use(condix)};
        nTrials = numel(trials2use);

        nEarly = sum(obj(sessix).bp.early(trials2use));
        early_per_trial(sessix,condix) = nEarly / nTrials;
    end
end

f = figure; 
f.Position = [316          73        1232         905];
ax = axes(f);
nCond = numel(cond2use);
violincols = reshape(cell2mat(dfparams.plt.color(cond2use)),3,nCond)';
vs = violinplot(early_per_trial,dfparams.cond(cond2use),...
    'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1}, 'ViolinColor', violincols);
ylabel('fraction of early licks per trial')
ylim([0,1])
title('early licks all sessions, all mice')
ax.FontSize = 20;

%% kinematics

[kin,kinfeats] = getKin(meta,obj,dfparams,params);

%% plot kinematics

% feats2plot = {'tongue_ydisp_view1',...
%               'tongue_yvel_view1',...
%               'jaw_ydisp_view1',...
%               'jaw_yvel_view1'};

% feats2plot = {'tongue_ydisp_view1',...
%     'jaw_ydisp_view1',...
%     'jaw_yvel_view1'};
feats2plot = {'motion_energy'};
cond2plot = 3:6;
sav = 0;

plotKinfeats(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)

%% avg jaw velocity during stim

if strcmpi(dfparams.alignEv,'delay') % function depends on data aligned to delay period, since we use the stim period to measure avgjawvel
%     feats2plot = {'jaw_ydisp_view1',...
%                   'jaw_yvel_view1'};
    feats2plot = {'motion_energy'};
    cond2plot = 3:6;
    sav = 0;

    plotAvgJawVelocityDuringStim(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)
end

%% decode trial type from single trial jaw pos/vel

k = 1; % number of iterations (bootstrap)

dt = 4; % train/test a model every dt'th time point
mdlTime = dfparams.time(1:dt:numel(dfparams.time));
numT = numel(mdlTime);


% feats2use = {'jaw_ydisp_view1',...
%              'jaw_yvel_view1'};
% feats2use = {'jaw_ydisp_view1'};
feats2use = {'motion_energy'};
[~,mask] = patternMatchCellArray(kin(1).featLeg,feats2use,'any');
featix = find(mask);

train = 0.8; % fraction of trials to use for training (1-train for testing)


cond2use = [3 5]; % right,left
acc_nostim = kinChoiceDecoder(meta,numT,k,kinfeats,cond2use,params,train,featix,dt);

cond2use = [4 6]; % right stim,left stim
acc_stim = kinChoiceDecoder(meta,numT,k,kinfeats,cond2use,params,train,featix,dt);

% figure; imagesc(squeeze(acc_nostim)'); colorbar;
% figure; imagesc(squeeze(acc_stim)'); colorbar;


%% plots


align = mode(obj(1).bp.ev.(dfparams.alignEv));
sample = mode(obj(1).bp.ev.sample) - align;
delay = mode(obj(1).bp.ev.delay) - align;

f = figure; ax = axes(f); hold on;
stderr = std(std(acc_nostim,[],2),[],3) ./ (k*numel(meta));
means = mean(mean(acc_nostim,2),3);
shadedErrorBar(mdlTime, means, stderr, {'Color',dfparams.plt.color{1},'LineWidth',1.5},1, ax);

stderr = std(std(acc_stim,[],2),[],3) ./ (k*numel(meta));
means = mean(mean(acc_stim,2),3);
shadedErrorBar(mdlTime, means, stderr, {'Color',dfparams.plt.color{2},'LineWidth',1.5},1, ax);

xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

title('Decode choice from jaw movements')
xlabel(['Time (s) from ' dfparams.alignEv])
ylabel('Accuracy')
ax = gca;
ax.FontSize = 20;
xlim([mdlTime(5) mdlTime(end)])
hold off

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig3/figs/choiceDecoder';
% fn = [meta.anm '_' meta.date '_' params.lfads_run '_' 'winsize_' num2str(sm)];
% mysavefig(f,pth,fn)


%% same as above but plotting each session individually

align = mode(obj(1).bp.ev.(dfparams.alignEv));
sample = mode(obj(1).bp.ev.sample) - align;
delay = mode(obj(1).bp.ev.delay) - align;

f = figure; ax = axes(f); hold on;
stderr = squeeze(std(acc_nostim,[],2)) ./ k;
means = squeeze(mean(acc_nostim,2));
for i = 1:size(means,2)
    shadedErrorBar(mdlTime, means(:,i), stderr(:,i), {'Color',dfparams.plt.color{1},'LineWidth',0.5},1, ax);
end

stderr = squeeze(std(acc_stim,[],2)) ./ k;
means = squeeze(mean(acc_stim,2));
for i = 1:size(means,2)
    shadedErrorBar(mdlTime, means(:,i), stderr(:,i), {'Color',dfparams.plt.color{3},'LineWidth',0.5},1, ax);
end

xline(sample,'k--','LineWidth',2)
xline(delay,'k--','LineWidth',2)
xline(0,'k--','LineWidth',2)

title('Decode choice from jaw movements')
xlabel(['Time (s) from ' dfparams.alignEv])
ylabel('Accuracy')
ax = gca;
ax.FontSize = 20;
xlim([mdlTime(5) mdlTime(end)])
hold off



