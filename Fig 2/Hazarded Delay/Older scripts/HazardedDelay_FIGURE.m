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
%params.condition(1) = {'R&hit&~stim.enable&~early'}; % right hits, no stim, aw off
%params.condition(end+1) = {'L&hit&~stim.enable&~early'}; % left hits, no stim, aw off

% Delay period length that you want to warp all delay lengths to
params.desiredDelay = 0.9000;       

% Different delay period lengths that were used in hazarded delay
params.delay(1) = 0.3000;
params.delay(2) = 0.6000;
params.delay(3) = 1.2000;
params.delay(4) = 1.8000;
% params.delay(5) = 2.4000;
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
disp('Loading Hazarded Delay behavior objs')
for i = 1:numel(meta)
    obj = objs{i};
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
%%  PLOT SUMMARY FIGURES FOR EACH SESSION
% %   Figure 1 = prob of jaw movement for each delay duration (averaged
% %   throughout the session)
% for gg = 1:numel(meta)
%     obj = objs{gg};
%     met = meta(gg);
% 
%     anm = met.anm;
%     date = met.date;
%     sesstitle = strcat(anm,date);  % Name/title for session
% 
%     delaylen = obj.bp.ev.goCue - obj.bp.ev.delay;       % Find the delay length for all trials
%     conditions = {1,2};
% 
%     met = getDelayTrialID(met,conditions,delaylen);     % Group the trials in each condition based on their delay length
% 
%     % Find the probability of jaw [Jaw] movement at all time points in the session for trials of
%     % specific conditions
%     jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'prob',params);    % (1 x conditions cell array)
%     % Each cell: (time x trials in that condition)
% 
%     % Find average jaw velocity for each condition (all delay lengths
%     % included)
%     jawvel.right = mean(jaw_by_cond{1},2,'omitnan');         
%     jawvel.left = mean(jaw_by_cond{2},2,'omitnan');
%     
%     % Plot probability of jaw movement 
%     figure();
%     colors = {[0 0 1],[1 0 0]};
%     plot(taxis,jawvel.right,'Color',colors{1},'LineWidth',2)
%     hold on;
%     plot(taxis,jawvel.left,'Color',colors{2},'LineWidth',2)
%     legend('Right','Left')
%     xlim([-1.4 0.3])
%     sgtitle(sesstitle)
% end

%% Plot summary figure for average across all animals

%%%% Find average prob of jaw movement for each animal (across each of its
% sessions) %%%%
conditions = {1};
[jaw_allAnm,nAnimals,uc,sessbyAnm] = findJawProb_Multi(objs,meta,conditions,taxis,params);                  % Hazarded delay animals
[jaw_ctrl, nAnimals_ctrl,sessbyAnm_ctrl] = findJawProb_Multi(ctrlobjs,ctrlmeta,conditions,taxis,params);    % Static delay animals 

%%%% Average across all animals %%%%
% Hazarded delay animals
temp.haz = [];
%temp.left = [];
%temp.selectivity = [];
for a = 1:nAnimals
    temp.haz = [temp.haz,jaw_allAnm.haz{a}];
    %temp.left = [temp.left,jaw_allAnm.left{a}];
    %temp.selectivity = [temp.selectivity, jaw_allAnm.right{a}-jaw_allAnm.left{a}];
end

AvgAll.haz = medfilt1(mean(temp.haz,2,'omitnan'),10);
% AvgAll.left = medfilt1(mean(temp.left,2,'omitnan'),10);
% AvgAll.selectivity = mySmooth(mean(temp.selectivity,2,'omitnan'),51);
% Static delay animals 
Ctrltemp.haz = [];
% Ctrltemp.left = [];
% Ctrltemp.selectivity = [];
for a = 1:nAnimals_ctrl
    Ctrltemp.haz = [Ctrltemp.haz,jaw_ctrl.haz{a}];
%     Ctrltemp.left = [Ctrltemp.left,jaw_ctrl.left{a}];
%     Ctrltemp.selectivity = [Ctrltemp.selectivity, jaw_ctrl.right{a}-jaw_ctrl.left{a}];
end
AvgCtrl.haz = medfilt1(mean(Ctrltemp.haz,2,'omitnan'),10);
% AvgCtrl.left = medfilt1(mean(Ctrltemp.left,2,'omitnan'),10);
% AvgCtrl.selectivity = mySmooth(mean(Ctrltemp.selectivity,2,'omitnan'),51);

%%%% Plot probability of jaw movement %%%% 
figure();
colors = {[0 0 0],[0.5 0.5 0.5]};
plot(taxis,100*AvgAll.haz,'Color',colors{1},'LineWidth',3)
hold on;
%plot(taxis,100*AvgAll.left,'Color',colors{2},'LineWidth',3)
plot(taxis,100*AvgCtrl.haz,'Color',colors{2},'LineWidth',2.5,'LineStyle','-.')
%plot(taxis,100*AvgCtrl.left,'Color',colors{4},'LineWidth',3)
xlim([-1 2])
xline(0,'LineStyle','--')
xline(-1.3,'LineStyle','--')
legend('Hazarded delay','Static delay','Location','best')
sgtitle('Avg across all animals')
xlabel('Time before delay (s)')
ylabel('Probability of jaw movement (%)')

% figure();
% AvgAll.selectivity = AvgAll.selectivity./max(AvgAll.selectivity);
% AvgCtrl.selectivity = AvgCtrl.selectivity./max(AvgAll.selectivity);
% colors = {[0 0 1],[1 0 0],[0 0 0.67],[0.67 0 0]};
% plot(taxis,100*AvgAll.selectivity,'Color','green','LineWidth',3)
% hold on;
% plot(taxis,100*AvgCtrl.selectivity,'Color','magenta','LineWidth',3)
% xlim([-1.35 0.3])
% xline(0,'LineStyle','--')
% xline(-1.3,'LineStyle','--')
% legend('Rand selectivity','Static selectivity','Location','best')
% sgtitle('Avg across all animals')
% xlabel('Time before delay (s)')
% ylabel('Selectivity (R - L) (%)')




