% Quantifying behavioral performance and cortical dependence in the
% Alternating Context Task
% -------------------------------------------------------------------------------------
clear,clc,close all

% add paths
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig3')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 1';
addpath(genpath(fullfile(figpth,'Utils')));
%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% set conditions to calculate behavioral performance for 
params.condition(1)     = {'(hit|miss|no)'};                             % all trials

params.condition(end+1) = {'L&~stim.enable&~autowater&~early'};             % right, no stim, 2AFC
params.condition(end+1) = {'L&stim.enable&~autowater&~early'};              % right, stim, 2AFC
params.condition(end+1) = {'R&~stim.enable&~autowater&~early'};             % left, no stim, 2AFC
params.condition(end+1) = {'R&stim.enable&~autowater&~early'};              % left, stim, 2AFC

params.condition(end+1) = {'L&~stim.enable&autowater&~early'};              % right, no stim, AW
params.condition(end+1) = {'L&stim.enable&autowater&~early'};               % right, stim, AW
params.condition(end+1) = {'R&~stim.enable&autowater&~early'};              % left, no stim, AW
params.condition(end+1) = {'R&stim.enable&autowater&~early'};               % left, stim, AW
%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

meta = [];
meta = loadMAH13_MCStim(meta,datapth);
meta = loadMAH14_MCStim(meta,datapth);

%% LOAD DATA

% ----------------------------------------------
% -- Behavioral and Video Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadBehavSessionData(meta,params);
%% 
conds2use = 2:9;
lickCutoff = 0.6;
trialCutoff = 40;
blah = [];
for sessix = 1:length(meta)
    lastTrial = obj(sessix).bp.Ntrials-trialCutoff;
    for c = 1:length(conds2use)
        cond = conds2use(c);
        condtrix = params(sessix).trialid{cond};
        condtrix = condtrix(condtrix<lastTrial);
        lickInit = false(length(condtrix),1);
        for t = 1:length(condtrix)
            currtrial = condtrix(t);
            allLicks = [obj(sessix).bp.ev.lickL{currtrial},obj(sessix).bp.ev.lickR{currtrial}];
            allLicks = sort(allLicks,'ascend');
            if ~isempty(allLicks)
                go = obj(sessix).bp.ev.goCue(currtrial);
                postGoLicks = allLicks((allLicks>go));
                if ~isempty(postGoLicks)
                    firstLick_aligned = postGoLicks(1)-go;
                    if firstLick_aligned<lickCutoff && obj(sessix).bp.hit(currtrial)
                        lickInit(t) = 1;                     % 0 = no correct lick before LickCutoff time; 1 = there is a correct lick before LickCutoff
                    end
                end
            end
        end
        perf(c) = sum(lickInit)/length(condtrix);
    end
    rez(sessix).trixWLicks = perf;
end
%% Concatenate across sessions
clearvars -except obj meta rez params lickCutoff 
perf_all = [];
for sessix = 1:length(meta)
    perf_all = [perf_all;rez(sessix).trixWLicks];
end
%% INCLUSION CRITERIA: Omit sessions where 2AFC stim didn't work or where performance was bad
sess2omit = false(length(meta),1);
LAFCstimCond = 1;
RAFCstimCond = 3;
EffectCutoff = 0.15;
perfCutoff = 0.55;
for sessix = 1:length(meta)
    temp = perf_all(sessix,:);
    if temp(RAFCstimCond)-temp(RAFCstimCond+1) < EffectCutoff || temp(LAFCstimCond)-temp(LAFCstimCond+1) < EffectCutoff
        sess2omit(sessix) = 1;
    end
    bp = obj(sessix).bp;
    Rperf = sum(bp.R&bp.hit&~bp.autowater)/sum(bp.R&~bp.autowater);
    Lperf = sum(bp.L&bp.hit&~bp.autowater)/sum(bp.L&~bp.autowater);
    if Rperf<perfCutoff || Lperf<perfCutoff
        sess2omit(sessix) = 1;
    end
end
perf_all(sess2omit,:) = [];
%% Plot
clearvars -except obj meta rez params perf_all lickCutoff sess2omit


colors = {[1 0 0],[1 0 0],[0 0 1],[0 0 1]};
%plotPerfV1(colors,perf_all,lickCutoff)

colors = {[1 0 0],[1 0.5 0.5],[0 0 1],[0.5 0.5 1]};
plotPerfV2(colors,perf_all,lickCutoff)


%%
function plotPerfV1(colors,perf_all,lickCutoff)
figure();
subplot(1,2,1)
conds2plot = 1:4;
for x = conds2plot
    xx = x*ones(size(perf_all,1),1);
    if x==1||x==3
        scatter(xx,perf_all(:,x),'filled','MarkerFaceColor',colors{x}); hold on;
    else 
        scatter(xx,perf_all(:,x),'MarkerEdgeColor',colors{x}); hold on;
    end
    for sessix = 1:size(perf_all,1)
        plot(1:2,perf_all(sessix,1:2),'Color','black')
        plot(3:4,perf_all(sessix,3:4),'Color','black')
    end
end
xlim([0 5])
xticks([1,2,3,4])
xticklabels({'L ctrl', 'L stim','R ctrl', 'R stim'})
ylabel(['Proportion of trials w/ lick within ' num2str(lickCutoff) ' (s) of goCue'])
title('2AFC')

subplot(1,2,2)
conds2plot = 5:8;
for x = 1:length(conds2plot)
    cond = conds2plot(x);
    xx = x*ones(size(perf_all,1),1);
    if x==1||x==3
        scatter(xx,perf_all(:,cond),'filled','MarkerFaceColor',colors{x}); hold on;
    else 
        scatter(xx,perf_all(:,cond),'MarkerEdgeColor',colors{x}); hold on;
    end
    for sessix = 1:size(perf_all,1)
        plot(1:2,perf_all(sessix,conds2plot(1:2)),'Color','black')
        plot(3:4,perf_all(sessix,conds2plot(3:4)),'Color','black')
    end
end
xlim([0 5])
xticks([1,2,3,4])
xticklabels({'L ctrl', 'L stim','R ctrl', 'R stim'})
ylabel(['Proportion of trials w/ lick within ' num2str(lickCutoff) ' (s) of waterDrop'])
title('Autowater')
end

function plotPerfV2(colors,perf_all,lickCutoff)
figure();
subplot(1,2,1)
conds2plot = 1:4;
means = mean(perf_all,1,'omitnan');
for x = conds2plot
    if x==1||x==2
        scatter(x,means(:,x),65,'filled','MarkerFaceColor',colors{1}); hold on;
    else
        scatter(x,means(:,x),65,'filled','MarkerFaceColor',colors{3}); hold on;
    end
    plot(1:2,means(1:2),'Color',colors{1},'LineWidth',3)
    plot(3:4,means(3:4),'Color',colors{3},'LineWidth',3)
    for sessix = 1:size(perf_all,1)
        plot(1:2,perf_all(sessix,1:2),'Color',colors{2})
        plot(3:4,perf_all(sessix,3:4),'Color',colors{4})
    end
end
xlim([0 5])
xticks([1,2,3,4])
xticklabels({'L ctrl', 'L stim','R ctrl', 'R stim'})
ylabel('Proportion of trials')
title('2AFC')

subplot(1,2,2)
conds2plot = 5:8;
for x = 1:length(conds2plot)
    cond = conds2plot(x);
    if x==1||x==2
        scatter(x,means(:,cond),65,'filled','MarkerFaceColor',colors{1}); hold on;
    else
        scatter(x,means(:,cond),65,'filled','MarkerFaceColor',colors{3}); hold on;
    end
    plot(1:2,means(5:6),'Color',colors{1},'LineWidth',3)
    plot(3:4,means(7:8),'Color',colors{3},'LineWidth',3)
    for sessix = 1:size(perf_all,1)
        plot(1:2,perf_all(sessix,5:6),'Color',colors{2})
        plot(3:4,perf_all(sessix,7:8),'Color',colors{4})
    end
end
xlim([0 5])
xticks([1,2,3,4])
xticklabels({'L ctrl', 'L stim','R ctrl', 'R stim'})
%ylabel(['Frac. of trials w/ lick within ' num2str(lickCutoff) ' (s) of waterDrop'])
title('Autowater')

sgtitle(['Proportion of trials w/ correct lick within ' num2str(lickCutoff) ' (s) of goCue/waterDrop'])
end

