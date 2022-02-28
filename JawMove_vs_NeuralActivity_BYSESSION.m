% SCRIPT FOR RELATING JAW MOVEMENTS TO NEURAL ACTIVITY ON A SESSION BY
% SESSION BASIS

% Will generate a figure for each session:
% Panel A: scatter plot of single trial choice mode projections vs. jaw
% velocity
% Panel B: Heatmap of jaw movements on each non-early, 2AFC hit trial
% Panel C: Projections of R and L hit, non-early, 2AFC trials onto choice
% modes found w/ and w/o early move trials
% Panel D: Avg Probability of jaw movement throughout the trial

% In Panels A and C--> choice mode is found according to Nuo's definition (hits
% and misses included)

%%
clear; clc; close all;

% Saving params
outputdir = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Figures\Uninstructed Movements';
toSave = 'yes';
%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % goCue or firstLick

params.lowFR               = 1; % remove clusters firing less than this val

% set conditions to use for projections
params.condition(1) = {'R&hit&~stim.enable&autowater.nums==2&~early'}; % right hits, no stim, aw off
params.condition(2) = {'L&hit&~stim.enable&autowater.nums==2&~early'}; % left hits, no stim, aw off
params.condition(3) = {'R&miss&~stim.enable&autowater.nums==2&~early'};   % error right, no stim, aw off
params.condition(4) = {'L&miss&~stim.enable&autowater.nums==2&~early'};   % error left, no stim, aw off
params.condition(5) = {'R&hit&~stim.enable&autowater.nums==1&~early'}; % right hits, no stim, aw on
params.condition(6) = {'L&hit&~stim.enable&autowater.nums==1&~early'}; % left hits, no stim, aw on
params.condition(7) = {'~hit&~miss&~stim.enable&autowater.nums==2&~early'}; % ignore, 2afc, no stim


% set conditions used for finding the modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};     % R hits, 2afc, stim on/off, not early
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};     % L hits, 2afc, stim on/off, not early
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};    % R miss, 2afc, stim on/off, not early
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};    % L miss, 2afc, stim on/off, not early
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};       % All hits, 2afc, stim on/off, not early
params.modecondition(6) = {['hit&autowater.nums==1&stim.num==' stim '&~early']};        % All hits, aw on, stim on/off, not early

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
% meta = loadJEB7_ALMVideo(meta);
meta = loadEKH1_ALMVideo(meta);
meta = loadEKH3_ALMVideo(meta);
meta = loadJGR2_ALMVideo(meta);
meta = loadJGR3_ALMVideo(meta);

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

for gg = 1:length(meta)
    f = figure(gg);
    f.WindowState = 'maximized';
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


    % PANEL A: SINGLE-TRIAL CHOICE MODE (rez) VS. JAW VELOCITY
    % This choice mode is found using R and L hits and misses 
    allModes = calcAllModes(obj,met,rez,params);
      
    % Project single trials onto choice mode
    cd = allModes.choice_mode;
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
    subplot(3,2,6)
    ActivityMode_Jaw_Scatter(jawVel_late,Choice_late,conditions,met,clrs)
    lgd = legend('Right','Left');
    lgd.FontSize = 12;

    % PANEL B: HEATMAP OF JAW MOVEMENTS THROUGHOUT SESSION
    conditions = {1,2};
    subplot(3,2,5)
    JawVelHeatmap(conditions,jaw,taxis,met)

    % PANEL C: CHOICE MODE WITH AND WITHOUT EARLY MOVEMENT TRIALS
    % This choice mode is only found using L and R hit trials, epoch =
    % delay
      
    % Find early tongue and jaw movements 
    obj = findEarlyMoveTrials(obj);

    % Find all modes when early trials are not included
    allModes_NoMove = calcAllModes(obj,met,rez,params);

    % Which conditions to project onto the modes
    conditions = [1,2];         % Left and right 2AFC hits (not early)
    smooth = 61;
    [allModes,allModes_NoMove] = getChoiceModeProjections_twoModes(obj,allModes,allModes_NoMove,smooth,conditions);

    % Plot the R and L hits onto the two choice modes (found w/ and w/o
    % early move trials) 
    if ~isempty(allModes.latentChoice) && ~isempty(allModes_NoMove.latentChoice)
        subplot(3,2,3)
        colors = {[0 0 1],[1 0 0],[0.5 0.5 1],[1 0.5 0.5]};
        lw = 2;
        plot(rez.time,allModes.latentChoice{1},'Color',colors{1},'LineWidth',2)
        hold on
        plot(rez.time,allModes.latentChoice{2},'Color',colors{2},'LineWidth',2)
        plot(rez.time,allModes_NoMove.latentChoice{1},'Color',colors{3},'LineWidth',2)
        plot(rez.time,allModes_NoMove.latentChoice{2},'Color',colors{4},'LineWidth',2)
        legend('Right','Left','R no early','L no early')
        title('Remove early move trials from choice mode','FontSize',14)
        xlabel('Time since go-cue (s)','FontSize',13)
        ylabel('Choice Mode (a.u.)','FontSize',13)
        xlim([-2.5 2.5])
    end
   

    % PANEL D: PLOT PROB OF JAW MOVEMENT ACROSS TRIAL FOR R AND L HITS
    subplot(3,2,2)
    dirs = {'R','L'};                           % Dimension of 'dirs' must match dimension of 'contexts'
    contexts = {'afc','afc'};           
    plotJawProb_SessAvg(dirs, contexts,anm,date,colors,edges)
    legend('Right','Left')

    sgtitle(sesstitle,'FontSize',16)

    % PANEL E: PROJECT NO-EARLYMOVE TRIALS ONTO REGULAR CHOICE MODE
    subplot(3,2,4)
    obj = removeEarlyTrials_PSTH(obj,met);              % Get trial PSTHs by condition, excluding early move trials
    [Reg_projection] = getChoiceModeProjections_oneMode(obj,allModes,smooth,conditions);      % Project all trials and non-early move trials onto same choice mode

    if ~isempty(Reg_projection.ALLTrix) && ~isempty(Reg_projection.NOEarly)
        colors = {[0 0 1],[1 0 0],[0 0 0.75],[0.75 0 0]};
        lw = 2;
        plot(rez.time,Reg_projection.ALLTrix{1},'Color',colors{1},'LineWidth',2)
        hold on
        plot(rez.time,Reg_projection.ALLTrix{2},'Color',colors{2},'LineWidth',2)
        plot(rez.time,Reg_projection.NOEarly{1},'Color',colors{3},'LineWidth',2)
        plot(rez.time,Reg_projection.NOEarly{2},'Color',colors{4},'LineWidth',2)
        legend('Right','Left','R no early','L no early')
        title('Remove early move trials from projection','FontSize',14)
        xlabel('Time since go-cue (s)','FontSize',13)
        ylabel('Choice Mode (a.u.)','FontSize',13)
        xlim([-2.5 2.5])
    end

    % PANEL F: Heat-map of R/L selectivity for all cells
    subplot(3,2,1)
    numCells = size(obj.psth,2);
    Rpsth = obj.psth(:,:,1);                        
    Lpsth = obj.psth(:,:,2);
    DirSelect = Rpsth-Lpsth;
    maxSelect =  max(abs(DirSelect),[],'all');      % Find max directional selectivity
    selectNorm = DirSelect./maxSelect;              % Normalize all selectivity values to max selectivity

    alignEv = find(taxis == 0);
    earlyDelay = find(taxis>-0.855 & taxis<-0.845);     % t = 0.85 s before alignEvent
    lateDelay = find(taxis>-0.253 & taxis<-0.247);      % t = 0.25 s before alignEvent
    sort_by = lateDelay;          % alignEv, early delay, late delay, presample

    plotSelectivityHeatmap(taxis,sort_by,numCells,selectNorm)
    ylabel('Neuron #','FontSize',18)
    xlabel('Time since go cue (s)','FontSize',18)
    ax = gca;
    ax.FontSize = 14;
    title('Direction selectivity (R - L)','FontSize',18)
    hold off;

    if strcmp(toSave,'yes')
        saveas(gcf,fullfile(outputdir,sesstitle),'jpeg')
        close all
    end
end