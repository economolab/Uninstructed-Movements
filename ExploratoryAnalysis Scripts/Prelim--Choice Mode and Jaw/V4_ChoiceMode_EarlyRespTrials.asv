% SCRIPT FOR RELATING CHOICE MODE ACTIVITY TO JAW VELOCITY IN A GRADED
% MANNER (HOW DOES CHOICE MODE CHANGE WHEN JAW VELOCITY CHANGES?)

% Will generate summary figure relating jaw velocity to choice coding
% direction

% PARAMS THAT NEED TO BE SET:
% Saving params: 'toSave' --> whether or not to save figures
% Run params: 'earlytrials' --> which method you want to use to identify
% early movement trials
% 'moveThresh' --> what percentage of the delay period you want to be used
% for identifying early move trials
% 'alignEvent' --> which behavioral event you want to align the PSTHs and
% modes to
%%
clear; clc; close all;

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\ActivityModes'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Utils'));

% addpath(genpath('C:\Code\ActivityModes'));
% addpath(genpath('C:\Code\Uninstructed Movements\Uninstructed-Movements\DataLoadingScripts'));
% addpath(genpath('C:\Code\Uninstructed-Movements'));
% addpath(genpath('C:\Code\Utils'));


% Saving params
outputdir = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Figures\Uninstructed Movements';
toSave = 'no';
%% SET RUN PARAMS

% Which method you want to use to identify early movement trials:
% 'motionEnergy' or 'DeepLabCut'
params.earlytrials         =  'motionEnergy';
params.moveThresh          = 0.15;      % What percentage of the delay period you want to use for identifying early move trials
params.alignEvent          = 'goCue';   % goCue or firstLick
params.lowFR               = 1; % remove clusters firing less than this val
params.dt = 0.05;

% set conditions to use for projections
params.condition(1) = {'R&hit&~stim.enable&autowater.nums==2&~early'}; % right hits, no stim, aw off
params.condition(2) = {'L&hit&~stim.enable&autowater.nums==2&~early'}; % left hits, no stim, aw off
params.condition(3) = {'R&miss&~stim.enable&autowater.nums==2&~early'};   % error right, no stim, aw off
params.condition(4) = {'L&miss&~stim.enable&autowater.nums==2&~early'};   % error left, no stim, aw off
params.condition(5) = {'R&hit&~stim.enable&autowater.nums==1&~early'}; % right hits, no stim, aw on
params.condition(6) = {'L&hit&~stim.enable&autowater.nums==1&~early'}; % left hits, no stim, aw on
params.condition(7) = {'~hit&~miss&~stim.enable&autowater.nums==2&~early'}; % ignore, 2afc, no stim
params.condition(8) = {'R&hit&~stim.enable&autowater.nums==2&early'}; % right EARLY RESPONSE hits, no stim, aw off
params.condition(9) = {'L&hit&~stim.enable&autowater.nums==2&early'}; % left EARLY RESPONSE hits, no stim, aw off


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
meta = loadJEB4_ALMVideo(meta);
meta = loadJEB5_ALMVideo(meta);
meta = loadJEB6_ALMVideo(meta);
meta = loadJEB7_ALMVideo(meta);
% meta = loadEKH1_ALMVideo(meta);
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
%%
for gg = 1:length(meta)         % For all loaded sessions...
    ff = figure(gg);
    %ff.WindowState = 'maximized';
    obj = objs{gg};
    met = meta(gg);

    anm = obj.pth.anm;                  % Animal name
    date = obj.pth.dt;                  % Session date
    probenum = string(met.probe);       % Which probe was used
    sesstitle = strcat(anm,date,' ;  ','Probe ',probenum,'LateDelay');  % Name/title for session

    clear rez; clear removeEarly, clear reg

    rez.time = objs{1}.time;
    rez.condition = objs{1}.condition;
    rez.alignEvent = params.alignEvent;

    % Find all modes
    allModes = calcAllModes(obj,met,rez,params,'no');

    % Project single trials onto choice mode
    cd = allModes.choice_mode;
    conditions = {1,2,8,9};
    latent = getTrialLatents(obj,cd,conditions,met);

    % Which conditions to project onto the modes
    conditions = [1,2,8,9];         % Left and right 2AFC hits (not early)
    smooth = 61;
    allModes = getChoiceModeProjection(obj,allModes,smooth,conditions);

    % PANEL A: Plot projections of trials from specific conditions onto
    % choice mode
    if ~isempty(allModes.latentChoice)
        colors = {[0 0 1],[1 0 0],[0.5 0.5 0.9],[0.9 0.5 0.5]};
        lw = 3;
        for i=1:length(conditions)
            plot(rez.time,allModes.latentChoice{i},'Color',colors{i},'LineWidth',lw)
            hold on;
        end
        hold off;

        legend('Right','Left','Right Early Resp','Left Early Resp','Location','best')
        title(sesstitle,'FontSize',14)
        xlabel('Time since go-cue (s)','FontSize',13)
        ylabel('Delay CD (a.u.)','FontSize',13)
        xlim([-2.5 2.5])
    end

    for i=1:length(conditions)
        cond = conditions(i);
        trialid = met.trialid{cond};
        for tt = 1:length(trialid)
            
        end
    end
end

%sgtitle(sesstitle,'FontSize',16)

if strcmp(toSave,'yes')
    saveas(gcf,fullfile(outputdir,sesstitle),'jpeg')
    close all
end