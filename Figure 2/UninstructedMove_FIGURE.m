%%
clear; clc; close all;
%%
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\ActivityModes'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Utils'));

% addpath(genpath('C:\Code\ActivityModes'));
% addpath(genpath('C:\Code\Uninstructed Movements\Uninstructed-Movements\DataLoadingScripts'));
% addpath(genpath('C:\Code\Uninstructed-Movements'));
% addpath(genpath('C:\Code\Utils'));
% addpath(genpath('C:\Code\DataLoadingScripts'));


% Saving params
outputdir = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Figures\Uninstructed Movements';
toSave = 'no';
%% SET RUN PARAMS

% Which method you want to use to identify early movement trials:
% 'motionEnergy' or 'DeepLabCut'
params.alignEvent          = 'goCue';   % goCue or firstLick
params.lowFR               = 1; % remove clusters firing less than this val
params.dt = 0.05;
params.jawMeasure          = 'sideJaw'; % sideJaw or Trident

% set conditions to use for projections
params.condition(1) = {'R&hit&~stim.enable&autowater.nums==2&~early'}; % right hits, no stim, aw off, no early response
params.condition(2) = {'L&hit&~stim.enable&autowater.nums==2&~early'}; % left hits, no stim, aw of, no early response
% params.condition(3) = {'R&miss&~stim.enable&autowater.nums==2&~early'};   % error right, no stim, aw off
% params.condition(4) = {'L&miss&~stim.enable&autowater.nums==2&~early'};   % error left, no stim, aw off
% params.condition(5) = {'R&hit&~stim.enable&autowater.nums==1&~early'}; % right hits, no stim, aw on
% params.condition(6) = {'L&hit&~stim.enable&autowater.nums==1&~early'}; % left hits, no stim, aw on
% params.condition(7) = {'~hit&~miss&~stim.enable&autowater.nums==2&~early'}; % ignore, 2afc, no stim
% params.condition(8) = {'R&hit&~stim.enable&autowater.nums==2&early'}; % right EARLY RESPONSE hits, no stim, aw off
% params.condition(9) = {'L&hit&~stim.enable&autowater.nums==2&early'}; % left EARLY RESPONSE hits, no stim, aw off


% set conditions used for finding the modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};     % R hits, 2afc, stim on/off, not early
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};     % L hits, 2afc, stim on/off, not early
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};    % R miss, 2afc, stim on/off, not early
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};    % L miss, 2afc, stim on/off, not early
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};       % All hits, 2afc, stim on/off, not early
params.modecondition(6) = {['hit&autowater.nums==1&stim.num==' stim '&~early']};        % All hits, aw on, stim on/off, not early

%% SET METADATA FROM ALL RELEVANT SESSIONS/ANIMALS
meta = [];
% meta = loadJEB4_ALMVideo(meta);
% meta = loadJEB5_ALMVideo(meta);
% meta = loadJEB6_ALMVideo(meta);
% meta = loadJEB7_ALMVideo(meta);
meta = loadEKH1_ALMVideo(meta);
% meta = loadEKH3_ALMVideo(meta);
% meta = loadJGR2_ALMVideo(meta);
% meta = loadJGR3_ALMVideo(meta);

taxis = meta(end).tmin:meta(end).dt:meta(end).tmax;   % get time-axis with 0 as time of event you aligned to
taxis = taxis(1:end-1);
%% PREPROCESS DATA
objs = loadObjs(meta);

for i = 1:numel(meta)
    obj = objs{i};
    obj.condition = params.condition;
    % get trials and clusters to use
    meta(i).trialid = findTrials(obj, obj.condition);   % Get which trials pertain to the behavioral conditions you are looking at
    cluQuality = {obj.clu{meta(i).probe}(:).quality}';  % Get clusters that are of the qualities that you specified
    meta(i).cluid = findClusters(cluQuality, meta(i).quality);
    % align data
    obj = alignSpikes(obj,meta(i),params);              % Align the spike times to the event that you specified
    % get trial avg psth, single trial data, and single trial data grouped
    % by condition (aka R 2AFC, R AW, etc.)
    obj = getPSTHs(obj,meta(i));
    objs{i} = obj;
end
%% EXAMPLE HEATMAP OF JAW VELOCITY ON SINGLE TRIALS--SEPARATED BY TRIAL TYPE

sesh = 4;
obj = objs{sesh};     
met = meta(sesh);

anm = obj.pth.anm;                  % Animal name
date = obj.pth.dt;                  % Session date
probenum = string(met.probe);       % Which probe was used


% Find the jaw velocity at all time points in the session for trials of
% specific conditions
conditions = {1,2};
if strcmp(params.jawMeasure,'sideJaw')
    jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'vel');
elseif strcmp(params.jawMeasure,'Trident')
    jaw_by_cond = findTridentVelocity(taxis, obj,conditions,met);
end
l1 = size(jaw_by_cond{1},2);      % Number of trials in the first condition

% Sort the trials by jaw velocity during the late delay period
% Find the average jaw velocity during specified time points (on each trial)
    startix = find(taxis>=-0.4, 1, 'first');
    stopix = find(taxis<=-0.05, 1, 'last');
    val.right = nanmean(jaw_by_cond{1}(startix:stopix, :), 1);
    val.left = nanmean(jaw_by_cond{2}(startix:stopix, :), 1);
    nanix = find(isnan(val.right)); nanix = find(isnan(val.left));                   % Get rid of trials where jaw velocity is always NaN
    val.right(nanix) = [];  val.left(nanix) = [];
    
    % Sort the average jaw velocities in descending order and save the trial
    % order
    sort_by_cond = cell(1,numel(conditions));
    [sort_by_cond{1}, six1] = sort(val.right, 'descend'); jaw_by_cond{1} = jaw_by_cond{1}(:,six1);
    [sort_by_cond{2}, six2] = sort(val.left, 'descend'); jaw_by_cond{2} = jaw_by_cond{2}(:,six2);

% Plot
figure();
numTrixPlot = 30;
rangetoPlot = 1:numTrixPlot;
tGo = find(taxis==0);
for i=1:numel(conditions)
    if i==1
        subplot(1,2,2);
    elseif i==2
        subplot(1,2,1);
    end
    imagesc(taxis,1:numTrixPlot,jaw_by_cond{i}(:,1:numTrixPlot)'); caxis([0 1.5]); colormap("hot");
    go = 0;
    delstart = -0.9;
    sampstart = delstart-1.3;
    line([go,go],[0,numTrixPlot+0.5],'Color','white','LineStyle','--')
    line([delstart,delstart],[0,numTrixPlot+0.5],'Color','white','LineStyle','--')
    line([sampstart,sampstart],[0,numTrixPlot+0.5],'Color','white','LineStyle','--')
    xlim([taxis(10), taxis(tGo)])
    xlabel('Time before go-cue (s)','FontSize',13)
    ylabel('Trials','FontSize',13)
    if i==1
        title('Right trials')
    elseif i==2
        title('Left trials')
    end
    c=colorbar;
    ylabel(c,'Jaw velocity','FontSize',12,'Rotation',90);
end
figtitle =  strcat('Example Trials from',anm,date,' ;  ','Probe ',probenum);  % Name/title for session
sgtitle(figtitle,'FontSize',16)
%% EXAMPLE SESSION OF PROBABILITY OF JAW MOVEMENT--SEPARATED BY TRIAL TYPE

obj = objs{sesh};     % 11th data object = JEB7, 04-29 (Classic sesh)
met = meta(sesh);

anm = obj.pth.anm;                  % Animal name
date = obj.pth.dt;                  % Session date
probenum = string(met.probe);       % Which probe was used

conditions = {1,2};
colors = {[0 0 1],[1 0 0]};

% Plot
figure();
if strcmp(params.jawMeasure,'sideJaw')
    plotJawProb_SessAvg(obj,met,conditions,colors,taxis,'no',params)
elseif strcmp(params.jawMeasure,'Trident')
    plotTridentVel_SessAvg(obj,met,conditions,colors,'no',params)
end
legend('Right','Left','Location','best')

% Add lines at trial times
go = 0;
delstart = -0.9;
sampstart = delstart-1.3;
xline(go,'Color','black','LineStyle','--')
xline(delstart,'Color','black','LineStyle','--')
xline(sampstart,'Color','black','LineStyle','--')
xlim([-2.3 0])

figtitle =  strcat('Example w/ jaw selectivity',anm,date,' ;  ','Probe ',probenum);  % Name/title for session
title(figtitle,'FontSize',16)
%% EXAMPLE SESSION OF CHOICE MODE PROJECTIONS, SEPARATED BY TRIAL TYPE
obj = objs{sesh};     % 11th data object = JEB7, 04-29 (Classic sesh)
met = meta(sesh);

anm = obj.pth.anm;                  % Animal name
date = obj.pth.dt;                  % Session date
probenum = string(met.probe);       % Which probe was used

rez.time = objs{1}.time;
rez.condition = objs{1}.condition;
rez.alignEvent = params.alignEvent;

% Find CDchoice (coding dimension during delay period)
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
epoch = 'latedelay';
choice_mode = choiceMode(obj,met,cond,epoch,rez.alignEvent,'no');

% Get the projection of specified conditions onto the choice mode
smooth = 51;
conditions = [1,2];
latentChoice = getChoiceModeProjection(obj,choice_mode,smooth,conditions);

% Plot
figure();
colors = {[0 0 1],[1 0 0]};
for c = 1:numel(conditions)
    temp = conditions(c);
    col = colors{c};
    plot(taxis,latentChoice{temp},'Color',col,'LineWidth',3)
    hold on;
end
col = 'black';
addTrialLines(col)
xlim([-2.3 0])
xlabel('Time since go-cue')
ylabel('a.u.')
set(gca, 'YDir','reverse')
legend('Right','Left','Location','best')

figtitle =  strcat('Example CDchoice',anm,date,' ;  ','Probe ',probenum);  % Name/title for session
title(figtitle,'FontSize',16)
%% EXAMPLE SESSION--SCATTER PLOT OF SINGLE TRIAL JAW VELOCITY VS CHOICE CD

obj = objs{17};     % 11th data object = JEB7, 04-29 (Classic sesh)
met = meta(17);

anm = obj.pth.anm;                  % Animal name
date = obj.pth.dt;                  % Session date
probenum = string(met.probe);       % Which probe was used

rez.time = objs{1}.time;
rez.condition = objs{1}.condition;
rez.alignEvent = params.alignEvent;

% Project single trials onto choice mode
cd = choice_mode;
latent = getTrialLatents(obj,cd,conditions,met);
lat_choice = [];
jaw = [];
for c = 1:numel(conditions)
    lat_choice = [lat_choice,latent{c}];
    jaw = [jaw,jaw_by_cond{c}];
end

% Define time intervals: Time frame for late delay period(from -0.4 before go-cue to -0.1)
late_start = find(taxis>=-0.4, 1, 'first');
late_stop = find(taxis<=-0.05, 1, 'last');
lateDelay = late_start:late_stop;

% Get jaw velocity and activity mode averages for late delay
timeInt = lateDelay;
jawVel_late = getAverages(timeInt,jaw);
Choice_late = getAverages(timeInt,lat_choice);

% Make scatter plot
conditions = 1:2;               % Look only at correct left and right hits during 2AFC
figure();
ActivityMode_Jaw_Scatter(jawVel_late,Choice_late,conditions,met,colors,obj,params);

nanix = find(isnan(jawVel_late));                    % Find indices where the jaw vel is a NaN
jv = jawVel_late(~isnan(jawVel_late));               % Get rid of the NaN values
Choice_late(nanix) = [];                             % Indices that were a NaN for jaw vel, get rid of those indices in the choice mode as well
ch = Choice_late;
R = corr2(jv,ch);                           % Calculate the correlation coefficient between these two variables
R = num2str(R);
coeff = polyfit(jv,ch,1);                   % Find the line of best fit
hline = refline(coeff);
hline.LineStyle = '--';
hline.Color = 'k';
str = strcat('R^2 =',R);
lgd = legend('Right','Left',str);
lgd.FontSize = 11;
lgd.Location = 'best';
xlabel('Avg Jaw Velocity','fontsize',14)
xlim([0 1.2])
ylabel('Avg choice mode','fontsize',14)
figtitle =  strcat('Example Single trial correlations',anm,date,' ;  ','Probe ',probenum);  % Name/title for session
sgtitle(figtitle)

disp('hi');
            
%% PROB OF JAW MOVEMENT--ACROSS ALL ANIMALS
% conditions = {1,2};
% 
% jaw_allAnm = jawProb_ALLAnimals(objs,meta,taxis, conditions);   % (animals x conditions) cell array.  Each cell contains the average prob of jaw movement for that animal across sessions for that condition
% 
% figure();
% for c = 1:numel(conditions)
%     col = colors{c};
%     meanJaw = NaN(length(taxis)+1,nAnimals);
%     for a = 1:nAnimals
%         meanJaw(:,a) = jaw_allAnm{a,c};
%     end
%     meanJaw = mean(meanJaw,2,'omitnan'); meanJaw = meanJaw(2:end);
%     plot(taxis,meanJaw,'Color',col,'LineWidth',3)
%     hold on;
% end
% col = 'black';
% addTrialLines(col)
% legend('Right','Left')
% xlim([-2.3 2.5])
% xlabel('Time since go-cue (s)')
% ylabel('Probability of jaw movement')
% title('Prob of jaw movement -- ALL ANIMALS')


% %% CHOICE MODE FOR ALL SESSIONS AND ANIMALS
% multipsth = concatPSTH(objs);
% 
% rez.time = objs{1}.time;
% rez.condition = objs{1}.condition;
% rez.alignEvent = params.alignEvent;
% 
% % Find CDchoice (coding dimension during delay period)
% cond{1} = params.modecondition{1};
% cond{2} = params.modecondition{2};
% epoch = 'delay';
% choice_mode = choiceModeMulti(objs,meta,cond,epoch,rez.alignEvent,'no');
% 
% % Get the projection of specified conditions onto the choice mode
% smooth = 51;
% conditions = [1,2];
% latentChoice = getChoiceModeProjection_Multi(multipsth,choice_mode,smooth,conditions);
% 
% figure();
% colors = {[0 0 1],[1 0 0]};
% for c = 1:numel(conditions)
%     temp = conditions(c);
%     col = colors{c};
%     plot(taxis,latentChoice{temp},'Color',col,'LineWidth',3)
%     hold on;
% end
% col = 'black';
% addTrialLines(col)
% xlim([-2.3 2.5])
% xlabel('Time since go-cue')
% ylabel('a.u.')
% set(gca, 'YDir','reverse')
% legend('Right','Left','Location','best')
% 
% title('CDchoice - ALL Animals','FontSize',16)