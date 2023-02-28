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
params.condition(1) = {'R&hit&~stim.enable&~early'};  % R hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~early'};

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
%%

%%%% Find average prob of jaw movement for each animal (across each of its
% sessions) %%%%
conditions = {1,2};

[animNames,uc,nAnimals] = getAnimalNames(meta);
sessbyAnm = cell(1,nAnimals);
nConditions = numel(conditions);

jaw_allAnmhaz = NaN(length(taxis),length(params.delay),nConditions, length(meta)); % time x delay lengths x conditions
cnt = 0;
for a = 1:nAnimals
    currAnm = uc(a);    % Get the current animal name
    temp = strcmp(animNames,currAnm);     % Find the entries in 'meta' and 'obj' that correspond to this animal
    sessID = find(temp);
    nSess = length(sessID); sessbyAnm{a} = nSess;           % Number of sessions for this animal
    for i = 1:nSess               % For each session for this animal...
        cnt = cnt+1;
        currSess = sessID(i);
        met = meta(currSess);
        obj = objs{currSess};
        met.alignEvent = params.alignEvent;

        delaylen = obj.bp.ev.goCue - obj.bp.ev.delay;       % Find the delay length for all trials

        met = getDelayTrialID(met,conditions,delaylen);     % Group the trials in each condition based on their delay length

        % Find the probability of jaw [trident] movement at all time points in the session for trials of
        % specific conditions
        jaw_by_cond = findJawVelocity(taxis, obj,conditions,'prob',met);    % (1 x conditions cell array)
        % Each cell: (time x trials in that condition)
    
        for p = 1:length(params.delay)                  % For each delay length...
            for c = 1:nConditions
                gix = find(met.del_trialid{c}==p);              % Get the trial IDs in the first condition that have the current delay length
                tempjaw = mean(jaw_by_cond{c}(:,gix),2,'omitnan');     % Find avg jaw velocity for first condition trials with that delay
                jaw_allAnmhaz(:,p,c,cnt) = tempjaw;
            end
        end
    end

end

taxis_haz = taxis+0.5;              % Have to adjust for the time-axis offset in the data without SpikeGLX
del2use = 3;
All.avg.R = medfilt1(mean(jaw_allAnmhaz(:,del2use,1,:),4,'omitnan'),10);
All.avg.L = medfilt1(mean(jaw_allAnmhaz(:,del2use,2,:),4,'omitnan'),10);
All.std.R = std(jaw_allAnmhaz(:,del2use,1,:),0,4,'omitnan');
All.std.L = std(jaw_allAnmhaz(:,del2use,2,:),0,4,'omitnan');


%% Plot summary figure for average across all animals
%%%% Plot probability of jaw movement %%%% 
figure();
colors = {[0 0 1],[1 0 0]};
alpha = 0.2;
for c = 1:nConditions
    ax = gca;
    if c==1
        fn = 'R';
    else
        fn = 'L';
    end
    toplot = All.avg.(fn);
    err = 100*(1.96*(All.std.(fn)/length(meta)));
    shadedErrorBar(taxis_haz, 100*toplot, err ,{'Color',colors{c},'LineWidth',2}, alpha, ax); hold on;
end
xlim([-1.3 1.25])
xline(0,'LineStyle','--','LineWidth',1.5,'Color',[0 0 0])
xline(1.2 ,'LineStyle','-.','LineWidth',1.5,'Color',[0 0 0])
xline(-1.2,'LineStyle','--','LineWidth',1.5,'Color',[0 0 0])
%legend('R','L','Location','best')
sgtitle('Probability of jaw move')
xlabel('Time from delay (s)')
ylabel('%')
