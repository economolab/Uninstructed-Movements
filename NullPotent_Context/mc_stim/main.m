clear, clc, close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
rmpath(genpath(fullfile(utilspth,'fig3/')))
rmpath(genpath(fullfile(utilspth,'fig1/')))


% add paths for figure specific functions
addpath(genpath(pwd))


%% TODO
% - created new objs from pipeline and nidq data - only use those
% - need to handle frameTimes appropriately now (subtract 0.5 secs)
% - for choice decoding, just use avg during stim period

%% default params (same for each session)

dfparams = [];

% -- time alignment params --
% dfparams.alignEv = 'goCue';
% dfparams.times = [-2.5 2.5]; % relative to goCue
dfparams.alignEv = 'delay';
dfparams.times = [-1.0 2.5]; % relative to delay

dfparams.dt_vid = 0.0025;
dfparams.time = dfparams.times(1):dfparams.dt_vid:dfparams.times(2);

dfparams.warp = 0; % 0 means no warping, 1 means warp delay period to 


% -- trial type params --
dfparams.cond(1) = {'(hit|miss|no)&~stim.enable&~autowater&~autolearn'}; % all trials, no stim, no autowater, no autolearn
dfparams.cond(end+1) = {'(hit|miss|no)&stim.enable&~autowater&~autolearn'};  % all trials trials, stim, no autowater, no autolearn
dfparams.cond(end+1) = {'R&~stim.enable&~autowater&~autolearn'}; % right trials, no stim, no autowater
dfparams.cond(end+1) = {'R&stim.enable&~autowater&~autolearn'};  % right trials, stim, no autowater
dfparams.cond(end+1) = {'L&~stim.enable&~autowater&~autolearn'}; % left trials, no stim, no autowater
dfparams.cond(end+1) = {'L&stim.enable&~autowater&~autolearn'};  % left trials, stim, no autowater
dfparams.cond(end+1) = {'R&hit&~stim.enable&~autowater&~autolearn'}; % right hit trials, no stim, no autowater
dfparams.cond(end+1) = {'R&hit&stim.enable&~autowater&~autolearn'};  % right hit trials, stim, no autowater
dfparams.cond(end+1) = {'L&hit&~stim.enable&~autowater&~autolearn'}; % left hit trials, no stim, no autowater
dfparams.cond(end+1) = {'L&hit&stim.enable&~autowater&~autolearn'};  % left hit trials, stim, no autowater
% dfparams.cond(end+1) = {'R&miss&~stim.enable&~autowater&~autolearn'}; % left hit trials, no stim, no autowater
% dfparams.cond(end+1) = {'L&miss&stim.enable&~autowater&~autolearn'};  % left hit trials, stim, no autowater

% -- stim types --
dfparams.stim.types = {'Bi_MC','Right_MC','Left_MC','Bi_ALM','Bi_M1TJ','Right_ALM','Right_M1TJ','Left_ALM','Left_M1TJ'}; % ALM_Bi is MC_Bi
% dfparams.stim.num   = logical([0 0 0 1 0 0 0 0 0]); % Bi_ALM
dfparams.stim.num   = logical([0 0 0 0 1 0 0 0 0]); % Bi_M1TJ

% -- plotting params --
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
meta = loadMAH14_MCStim(meta,datapth);


% subset based on stim types
stim2use = dfparams.stim.types(dfparams.stim.num);
use = false(size(meta));
for sessix = 1:numel(use)
    [~,mask] = patternMatchCellArray({meta(sessix).stimLoc}, stim2use,'any');
    if mask
        use(sessix) = true;
    end
end
meta = meta(use);


obj = loadObjs(meta);


%% find trials for each condition

for i = 1:numel(meta)
    params(i).trialid = findTrials(obj(i), dfparams.cond);
end

%%
if dfparams.warp
    obj = warpData(obj,params);
end

%% behavioral performance
close all

rez = getPerformance(meta,obj,params);

% plots
cond2use = 1:6;
connectConds = 0;
plotPerformanceAllMice(meta,obj,rez,dfparams,params,cond2use,connectConds)
% plotPerformanceEachMouse(meta,obj,rez,dfparams,params) % TODO


%% number of early licks per trial

% % TODO: classify early licks by kinematics, not bpod data
% 
% cond2use = 3:6;
% early_per_trial = nan(numel(obj),numel(cond2use)); % (sessions,conds)
% for sessix = 1:numel(obj)
%     for condix = 1:numel(cond2use)
%         trials2use = params(sessix).trialid{cond2use(condix)};
%         nTrials = numel(trials2use);
% 
%         nEarly = sum(obj(sessix).bp.early(trials2use));
%         early_per_trial(sessix,condix) = nEarly / nTrials;
%     end
% end
% 
% f = figure; 
% f.Position = [316          73        1232         905];
% ax = axes(f);
% nCond = numel(cond2use);
% violincols = reshape(cell2mat(dfparams.plt.color(cond2use)),3,nCond)';
% vs = violinplot(early_per_trial,dfparams.cond(cond2use),...
%     'EdgeColor',[1 1 1], 'ViolinAlpha',{0.2,1}, 'ViolinColor', violincols);
% ylabel('fraction of early licks per trial')
% ylim([0,1])
% title('early licks all sessions, all mice')
% ax.FontSize = 20;

%% kinematics

for sessix = 1:numel(obj)
    if ~isstruct(obj(sessix).me)
        temp = obj(sessix).me;
        obj(sessix).me = [];
        obj(sessix).me.data = temp;
        clear temp
    end
end
% kin.dat (don't use)
% kin.featLeg corresponds to 3rd dimension of kinfeats/kinfeats_norm
[kin,kinfeats,kinfeats_norm] = getKin(meta,obj,dfparams,params);

%% plot kinematics

% feats2plot = {'tongue_ydisp_view1',...
%               'tongue_yvel_view1',...
%               'jaw_ydisp_view1',...
%               'jaw_yvel_view1'};

% feats2plot = {'tongue_ydisp_view1',...
%     'jaw_ydisp_view2',...
%     'jaw_yvel_view2'};
% feats2plot = {'tongue_ydisp_view1','jaw_ydisp_view1','motion_energy'};
feats2plot = {'motion_energy'};
cond2plot = 1:2;
% cond2plot = 3:6;
% cond2plot = 7:10;
sav = 0;

plotKinfeats(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)



%% avg jaw velocity during stim
close all
if strcmpi(dfparams.alignEv,'delay') % function depends on data aligned to delay period, since we use the stim period to measure avgjawvel
%     feats2plot = {'jaw_ydisp_view1',...
%         'jaw_yvel_view1',...
%         'motion_energy'};
%     feats2plot = {'jaw_yvel_view1',...
%         'tongue_ydisp_view1'};
    feats2plot = {'motion_energy'};
    cond2plot = 3:6;
    sav = 0;

%     plotAvgJawVelocityDuringStim(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)
    plotAvgJawVelocityDuringStim_v2(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)
%     plotAvgJawVelocityDuringStim_firstHalfDelay(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)
%     plotAvgJawVelocityDuringStim_singleTrials(meta,obj,dfparams,params,kin,kinfeats,feats2plot,cond2plot,sav)

end

%% decode trial type from single trial jaw pos/vel

k = 3; % number of iterations (bootstrap)

binSize = 20; % ms
binSize = floor(binSize / (dfparams.dt_vid*1000)); % samples

mdlTime = dfparams.time(1:binSize:numel(dfparams.time));
numT = numel(mdlTime);


feats2use = {'jaw_ydisp_view1',...
             'jaw_yvel_view1',...
             'motion_energy'};
% feats2use = {'jaw_ydisp_view1'};
% feats2use = {'motion_energy'};
[~,mask] = patternMatchCellArray(kin(1).featLeg,feats2use,'any');
featix = find(mask);

train = 0.8; % fraction of trials to use for training (1-train for testing)


cond2use = [7 9]; % right hit,left hit,
acc_nostim = kinChoiceDecoder(meta,numT,k,kinfeats,cond2use,params,train,featix,binSize); % (time,bootiters,sessions)

cond2use = [8 10]; % right stim hit,left stim hit
acc_stim = kinChoiceDecoder(meta,numT,k,kinfeats,cond2use,params,train,featix,binSize);

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
shadedErrorBar(mdlTime, means, stderr, {'Color',dfparams.plt.color{3},'LineWidth',1.5},1, ax);

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



%% plot single trial traces for a feature (stim vs non-stim)
meix = find(ismember(kin(1).featLeg,'motion_energy'));
me = kinfeats{1}(:,:,meix);
trix1 = params(1).trialid{1};
trix2 = params(1).trialid{2};
nTrials = numel(trix2);
trix1 = randsample(trix1,nTrials);
[~,ix1] = min(abs(dfparams.time - 0));
[~,ix2] = min(abs(dfparams.time - 0.8));
figure; hold on
for i = 1:numel(trix1)
    subplot(1,2,1)
    hold on
    plot(dfparams.time(ix1:ix2)   ,me(ix1:ix2,trix1(i)) + i, 'k')
end
xlim([dfparams.time(ix1) dfparams.time(ix2)])
xlabel('Time (s) from delay')
ylabel('jaw_ydisp_view1','Interpreter','none')
ax = gca;
ax.FontSize = 15;

for i = 1:numel(trix2)
    subplot(1,2,2)
    hold on
    plot(dfparams.time(ix1:ix2)   ,me(ix1:ix2,trix2(i)) + i, 'b')
end
xlim([dfparams.time(ix1) dfparams.time(ix2)])
xlabel('Time (s) from delay')
ylabel('jaw_ydisp_view1','Interpreter','none')
ax = gca;
ax.FontSize = 15;