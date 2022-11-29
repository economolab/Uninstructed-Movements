% Script for quantifying jaw movements on hazarded delay 2AFC
%%
clear; clc; close all;

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\ActivityModes\funcs'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\functions'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\utils'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Utils'));

% Saving params
outputdir = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Figures\Uninstructed Movements';
toSave = 'no';
%% SET RUN PARAMS

% Which method you want to use to identify early movement trials:
% 'motionEnergy' or 'DeepLabCut'
params.alignEvent          = 'delay';   % goCue or firstLick
params.dt = 0.05;

% set conditions to use for projections
params.condition(1) = {'hit&~stim.enable&~early'};  % All hits, no stim, aw off

% Delay period length that you want to warp all delay lengths to
params.desiredDelay = 0.9000;       

% Different delay period lengths that were used in hazarded delay
params.delay(1) = 0.3000;
params.delay(2) = 0.6000;
params.delay(3) = 1.2000;
params.delay(4) = 1.8000;
params.delay(5) = 2.4000;
%%  LOAD META DATA
meta = [];
meta = loadJEB11_BehavVid(meta);
meta = loadJEB12_BehavVid(meta);
meta(end).tmin = -2.5; % (s) relative to params.alignEvent
meta(end).tmax = 3;  % (s) relative to params.alignEvent
meta(end).dt = 0.005;

taxis = meta(end).tmin:meta(end).dt:meta(end).tmax;   % get time-axis with 0 as time of event you aligned to
taxis = taxis(1:end-1);
%% PREPROCESS DATA
objs = loadBehavVid(meta);
for i = 1:numel(meta)
    obj = objs{i};
    disp(['Loading Hazarded Delay behavior obj for session ' num2str(i) ' out of ' num2str(numel(meta))])
    meta(i).trialid = findTrials(obj, params.condition);   % Get which trials pertain to the behavioral conditions you are looking at
end
%% SET METADATA FROM ALL RELEVANT SESSIONS/ANIMALS
ctrlmeta = [];
ctrlmeta = loadJEB6_ALMVideo(ctrlmeta);
ctrlmeta = loadJEB7_ALMVideo(ctrlmeta);
ctrlmeta = loadEKH1_ALMVideo(ctrlmeta);
ctrlmeta = loadEKH3_ALMVideo(ctrlmeta);
ctrlmeta = loadJGR2_ALMVideo(ctrlmeta);
ctrlmeta = loadJGR3_ALMVideo(ctrlmeta);

taxis = ctrlmeta(end).tmin:ctrlmeta(end).dt:ctrlmeta(end).tmax;   % get time-axis with 0 as time of event you aligned to
taxis = taxis(1:end-1);
%% PREPROCESS DATA
ctrlobjs = loadObjs(ctrlmeta);

for i = 1:numel(ctrlmeta)
    disp(['Pre-processing data for Static Delay session ' num2str(i) ' out of ' num2str(numel(ctrlmeta))])
    obj = ctrlobjs{i};
    obj.condition = params.condition;
    % get trials and clusters to use
    ctrlmeta(i).trialid = findTrials(obj, obj.condition);   % Get which trials pertain to the behavioral conditions you are looking at
    cluQuality = {obj.clu{ctrlmeta(i).probe}(:).quality}';  % Get clusters that are of the qualities that you specified
    ctrlmeta(i).cluid = findClusters(cluQuality, ctrlmeta(i).quality);
    % align data
    obj = alignSpikes(obj,ctrlmeta(i),params);              % Align the spike times to the event that you specified
    % get trial avg psth, single trial data, and single trial data grouped
    % by condition (aka R 2AFC, R AW, etc.)
    obj = getPSTHs(obj,ctrlmeta(i));
    ctrlobjs{i} = obj;
end
%% Plot summary figure for average across all animals

%%%% Find average prob of jaw movement for each animal (across each of its
% sessions) %%%%
conditions = {1};
[jaw_allAnm,nAnimals,uc,sessbyAnm] = findHazJaw_SingleCond_Multi(objs,meta,conditions,taxis,params);                  % Hazarded delay animals
[jaw_ctrl, nAnimals_ctrl,~,sessbyAnm_ctrl] = findJawProb_Multi(ctrlobjs,ctrlmeta,conditions,taxis,params);    % Static delay animals 

taxis_haz = taxis+0.5;              % Have to adjust for the time-axis offset in the data without SpikeGLX

%%%% Average across all animals %%%%
% Hazarded delay animals
temp.haz = [];
delaylen_toplot = 3;
for a = 1:nAnimals
    temp.haz = [temp.haz,jaw_allAnm.haz{a,delaylen_toplot}];
end
All.avg = medfilt1(mean(temp.haz,2,'omitnan'),10);
All.std = std(temp.haz,0,2,'omitnan');

% Static delay animals 
Ctrltemp.haz = [];
for a = 1:nAnimals_ctrl
    Ctrltemp.haz = [Ctrltemp.haz,jaw_ctrl.haz{a}];
end
Ctrl.avg = medfilt1(mean(Ctrltemp.haz,2,'omitnan'),10);
Ctrl.std = std(Ctrltemp.haz,0,2,'omitnan');

[upperci, lowerci] = getConfInt(meta, ctrlmeta,All, Ctrl);


%%
%%%% Plot probability of jaw movement %%%% 
figure();
colors = {[0 0 0],[0.5 0.5 0.5]};
plot(taxis_haz,100*All.avg,'Color',colors{1},'LineWidth',2);hold on;
plot(taxis,100*Ctrl.avg,'Color',colors{2},'LineWidth',1.5,'LineStyle','-.')
patch([taxis_haz(10:end) fliplr(taxis_haz(10:end))],[lowerci.haz(9:end)' fliplr(upperci.haz(9:end)')],[0.25 0.25, 0.25],'FaceAlpha',0.2,'EdgeColor','none')
patch([taxis(10:end) fliplr(taxis(10:end))],[lowerci.static(9:end)' fliplr(upperci.static(9:end)')],[0.5 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')
xlim([-1 1.2])
xline(0,'LineStyle','--','LineWidth',1.15,'Color',[0 0 0])
xline(0.9,'LineStyle','-.','LineWidth',1,'Color',[0.25 0.25 0.25])
%xline(1.2,'LineStyle','--','LineWidth',1,'Color',[0 0 0])
legend('Hazarded delay','Static delay','Location','best')
sgtitle('Avg across all animals')
xlabel('Time from delay (s)')
ylabel('Probability of jaw movement (%)')

%%
function [upperci, lowerci] = getConfInt(meta, ctrlmeta,All, Ctrl)
nSessions = length(meta);

upperci.haz = 100*(All.avg(2:end)+1.96*(All.std(2:end)/nSessions));  % Find the upper 95% confidence interval for each condition
lowerci.haz = 100*(All.avg(2:end)-1.96*(All.std(2:end)/nSessions));  % Find lower 95% condifence interval for each condition  
upperci.haz = fillmissing(upperci.haz,'next'); lowerci.haz = fillmissing(lowerci.haz,'next');

nSessions = length(ctrlmeta);
upperci.static = 100*(Ctrl.avg(2:end)+1.96*(Ctrl.std(2:end)/nSessions));  
lowerci.static = 100*(Ctrl.avg(2:end)-1.96*(Ctrl.std(2:end)/nSessions));  
upperci.static = fillmissing(upperci.static,'next'); lowerci.static = fillmissing(lowerci.static,'next');
end