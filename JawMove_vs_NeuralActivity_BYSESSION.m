% SCRIPT FOR RELATING JAW MOVEMENTS TO NEURAL ACTIVITY ON A SESSION BY
% SESSION BASIS

% Will generate a figure for each session:
% Panel A: scatter plot of single trial choice mode projections vs. jaw
% velocity
% Panel B: Heatmap of jaw movements on each non-early, 2AFC hit trial
% Panel C: Projections of R and L hit, non-early, 2AFC trials onto choice
% modes found w/ and w/o early move trials
% Panel D: Avg Probability of jaw movement throughout the trial

% In Panel A --> choice mode is found according to Nuo's definition (hits
% and misses included)
% In Panel C --> choice mode is found only using hits and misses

%%
clear; clc; close all;
%% SET RUN PARAMS
params.alignEvent          = 'goCue';

params.lowFR               = 1; % remove clusters firing less than this val

% set conditions to use for projections
params.condition(1) = {'R&hit&~stim.enable&autowater.nums==2&~early'}; % right hits, no stim, aw off
params.condition(2) = {'L&hit&~stim.enable&autowater.nums==2&~early'}; % left hits, no stim, aw off
params.condition(3) = {'R&miss&~stim.enable&autowater.nums==2&~early'};   % error right, no stim, aw off
params.condition(4) = {'L&miss&~stim.enable&autowater.nums==2&~early'};   % error left, no stim, aw off
% params.condition(5) = {'R&hit&~stim.enable&autowater.nums==1&~early'}; % right hits, no stim, aw on
% params.condition(6) = {'L&hit&~stim.enable&autowater.nums==1&~early'}; % left hits, no stim, aw on
% params.condition(7) = {'~hit&~miss&~stim.enable&autowater.nums==2&~early'}; % ignore, 2afc, no stim
% params.condition(8) = {'R&hit&~stim.enable&autowater.nums==2&early'}; % right early hits, no stim, aw off
% params.condition(9) = {'L&hit&~stim.enable&autowater.nums==2&early'}; % left early hits, no stim, aw off


% set conditions used for finding the modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};

clrs = cell(1,4);
clrs{1} = [0 0 1];  % Blue = hit, 2AFC, R
clrs{2} = [1 0 0];  % Red = hit, 2AFC, L
clrs{3} = [0.5 0.5 1];    % Light blue = (for R choice mode w/ early trials removed)
clrs{4} = [1 0.5 0.5];    % Light red = (for L choice mode w/ early trials removed)

%% SET METADATA FROM ALL RELEVANT SESSIONS/ANIMALS
meta = [];
% meta = loadJEB4_ALMVideo(meta);
% meta = loadJEB5_ALMVideo(meta);
% meta = loadJEB6_ALMVideo(meta);
meta = loadJEB7_ALMVideo(meta);
% meta = loadEKH1_ALMVideo(meta);
% meta = loadEKH3_ALMVideo(meta);
% meta = loadJGR2_ALMVideo(meta);
% meta = loadJGR3_ALMVideo(meta);

taxis = meta(end).tmin:meta(end).dt:meta(end).tmax;   % get time-axis with 0 as time of event you aligned to
taxis = taxis(1:end-1);
%% PREPROCESS DATA
objs = loadObjs(meta);

%% PREPROCESS DATA
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

[multi_psth,fromsess] = concatPSTH(objs);                          % Concatenate the trial-averaged PSTHs from all sessions into one structure 
%% ACTIVITY MODES
% After each mode is calculated according to its definition, the weights
% for each neuron (aka the coding direction for that mode) is stored in
% rez.mode

% To project single trials or trial-averaged PSTHs onto that mode, have to
% multiply PSTH by the coding direction

for gg = 14:length(meta)
    figure(gg)
    obj = objs{gg};
    met = meta(gg);

    anm = obj.pth.anm;                  % Animal name
    date = obj.pth.dt;                  % Session date
    probenum = string(met.probe);       % Which probe 
    sesstitle = strcat(anm,{'   '},date,' ;  ','Probe ',{'   '},probenum);  % Name for session

    clear rez; clear removeEarly, clear reg
    
    rez.time = objs{1}.time;
    rez.condition = objs{1}.condition;
    rez.alignEvent = params.alignEvent;

    cond{1} = params.modecondition{1};
    cond{2} = params.modecondition{2};
    cond{3} = params.modecondition{3};
    cond{4} = params.modecondition{4};
    epoch = 'delay';

    % PANEL A: SINGLE-TRIAL CHOICE MODE (rez) VS. JAW VELOCITY
    % This choice mode is found using R and L hits and misses 
    rez.choice_mode = choiceMode(obj,met,cond,epoch,rez.alignEvent);
    clear cond
    
    % Project single trials onto coding dimensions
    cd = rez.choice_mode;
    lat_choice = getTrialLatents(obj,cd);
    
    % Find the jaw velocity at all time points in the trial
    edges = rez.time;
    jaw = findJawVelocity(edges, obj);

    % Define time intervals: Time frame for late delay period(from -0.4 before go-cue to -0.1)
    late_start = find(taxis<-0.3995 & taxis>-0.4005);
    late_stop = find(taxis<-0.0995 & taxis>-0.1005);
    lateDelay = late_start:late_stop;

    % Get jaw velocity and activity mode averages for late delay
    timeInt = lateDelay;
    jawVel_late = getAverages(timeInt,jaw);
    Choice_late = getAverages(timeInt,lat_choice);

    % Make scatter plot
    conditions = 1:2;               % Look only at correct left and right hits during 2AFC
    subplot(2,2,1)
    ActivityMode_Jaw_Scatter(jawVel_late,Choice_late,conditions,met,clrs)
    lgd = legend('Right','Left');
    lgd.FontSize = 12;

    % PANEL B: HEATMAP OF JAW MOVEMENTS THROUGHOUT SESSION
    conditions = {1,2};
    JawVelHeatmap(conditions,jaw,taxis,met)

    % PANEL C: CHOICE MODE WITH (reg) AND WITHOUT EARLY MOVEMENT TRIALS
    % This choice mode is only found using L and R hit trials, epoch =
    % delay
    reg.time = obj.time;
    reg.psth = obj.psth;
    reg.condition = obj.condition;
    reg.alignEvent = params.alignEvent;

    %choice mode
    cond{1} = params.modecondition{1};
    cond{2} = params.modecondition{2};
    epoch = 'delay';
    reg.choice_mode = choiceMode_noError(obj,met,cond,epoch,rez.alignEvent);
   
    % Find early tongue and jaw movements 
    obj = findEarlyMoveTrials(obj);

    % Choice modes without early movements (removeEarly)
    removeEarly.time = obj.time;
    removeEarly.psth = obj.psth;
    removeEarly.condition = obj.condition;
    removeEarly.alignEvent = params.alignEvent;

    cond{1} = params.modecondition{1};
    cond{2} = params.modecondition{2};
    epoch = 'delay';
    removeEarly.choice_mode = choiceMode_noError(obj,met,cond,epoch,removeEarly.alignEvent);
    clear cond

    % Which conditions to project onto the modes
    conditions = [1,2];         % Left and right 2AFC hits (not early)
    smooth = 61;
    [rez,removeEarly] = getChoiceModeProjections(obj,rez,removeEarly,smooth,conditions);

    % Plot the R and L hits onto choice modes found w/ and w/o early move
    % trials 
    if ~isempty(rez.latentChoice) && ~isempty(removeEarly.latentChoice)
        subplot(2,2,3)
        colors = {[0 0 1],[1 0 0],[0.5 0.5 1],[1 0.5 0.5]};
        lw = 2;
        plot(rez.time,rez.latentChoice{1},'Color',colors{1},'LineWidth',2)
        hold on
        plot(rez.time,rez.latentChoice{2},'Color',colors{2},'LineWidth',2)
        plot(rez.time,removeEarly.latentChoice{1},'Color',colors{3},'LineWidth',2)
        plot(rez.time,removeEarly.latentChoice{2},'Color',colors{4},'LineWidth',2)
        legend('Right Hit','Left Hit','R no early','L no early')
        title('Remove early move trials','FontSize',14)
        xlabel('Time since go-cue (s)','FontSize',13)
        ylabel('Choice Mode (a.u.)','FontSize',13)
    end
   

    % PANEL D: PLOT PROB OF JAW MOVEMENT ACROSS TRIAL FOR R AND L HITS
    subplot(2,2,4)
    dirs = {'R','L'};                           % Dimension of 'dirs' must match dimension of 'contexts'
    contexts = {'afc','afc'};           
    plotJawProb_SessAvg(dirs, contexts,anm,date,colors)
    legend('Right','Left')

    sgtitle(sesstitle,'FontSize',16)


end

%% ANALYSIS FUNCTIONS



