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
% addpath(genpath('C:\Code\DataLoadingScripts'));


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
meta = loadEKH1_ALMVideo(meta);
meta = loadEKH3_ALMVideo(meta);
meta = loadJGR2_ALMVideo(meta);
meta = loadJGR3_ALMVideo(meta);

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
%% Remove unwanted sessions

% remove sessions with less than 40 trials of rhit and lhit each (same as
% hidehiko ppn paper)
use = false(size(objs));
for i = 1:numel(use)
    met = meta(i);
    check1 = numel(met.trialid{1}) > 40;
    check2 = numel(met.trialid{2}) > 40;
    if check1 && check2
        use(i) = true;
    end
end

meta = meta(use);
objs = objs(use);
%%
for gg = 1:length(meta)         % For all loaded sessions...
    ff = figure(gg);
    ff.WindowState = 'maximized';
    obj = objs{gg};
    met = meta(gg);
    
    anm = obj.pth.anm;                  % Animal name
    date = obj.pth.dt;                  % Session date
    probenum = string(met.probe);       % Which probe was used
       
    
    % FIND KINEMATIC MODES
    
    % get kinematics
    
    jawAngle = getJawAngle(taxis, obj, met);    
    kin = struct();
    
%     [kin.featPos,kin.featVel] = getFeatureKinematics(taxis,obj,conditions,met,view,feat);
    
    view = 1; % side
    feat = 2; % jaw
    [kin.jawPos,~] = getFeatureKinematics(taxis,obj,conditions,met,view,feat); 
    view = 1; % side
    feat = 3; % nose
    [kin.nosePos,~] = getFeatureKinematics(taxis,obj,conditions,met,view,feat);
    view = 1; % side
    feat = 1; % tongue
    [kin.tonguePos,~] = getFeatureKinematics(taxis,obj,conditions,met,view,feat);
    view = 2; % bottom
    feat = 5; % top_paw
    [kin.topPawPos,~] = getFeatureKinematics(taxis,obj,conditions,met,view,feat);
    view = 2; % bottom
    feat = 6; % bottom_paw
    [kin.bottomPawPos,~] = getFeatureKinematics(taxis,obj,conditions,met,view,feat);
    
    params.tix = 1:1100;       % time points to use when finding mode
    params.fcut = 50;          % smoothing cutoff frequency
    params.cond = 1:2;         % which conditions to use to find mode
    params.method = 'xcorr';   % 'xcorr' or 'regress' (basically the same)
    params.fa = false;          % if true, reduces neural dimensions to 10 with factor analysis
    
    % get modes based on single trial full neural data (or latents) and kinematic features 
    
    kinfns = fieldnames(kin);
    for i = 1:numel(kinfns)
        Y = kin.(kinfns{i}); % feature data to use to calculate mode
        [mode.(kinfns{i}), dat.(kinfns{i})] = findMode(obj, Y, params);
        proj.(kinfns{i}) = getProjection(dat.(kinfns{i}), mode.(kinfns{i}));
    end
  
    
    % orthogonalize modes
    
    kinmodes = zeros(numel(mode.(kinfns{1})),numel(kinfns));
    for i = 1:numel(kinfns)
        kinmodes(:,i) = mode.(kinfns{i});
    end
    % TODO: kinmodes has rank of 1, so orthModes has 1 column vector in it
    orthModes = gschmidt(kinmodes);
    
    for i = 1:numel(kinfns)
        orthmode.(kinfns{i}) = orthModes(:,i);
    end
    
    % project data onto orthmodes
    
    for i = 1:numel(kinfns)
        orthproj.(kinfns{i}) = getProjection(dat.(kinfns{i}), orthmode.(kinfns{i}));
    end
    
    % mode is a struct with unorthgonalized modes
    % orthmode is a struct with orthogonalized modes
    % proj is a struct with single trials projected onto modes
    % orthproj is a struct with single trials projected onto orthmodes
    
    % plot everything
    for i = 1:numel(kinfns)
        f(i) = figure(i); 
        imagesc(orthproj.(kinfns{i})')
        title(kinfns{i});
    end
    
    % next steps:
    % Qpotent = orthModes;
    % Qnull is then the PCs of the residual activity after removing Qpotent
    % from single trials (or PSTHs)    
end
