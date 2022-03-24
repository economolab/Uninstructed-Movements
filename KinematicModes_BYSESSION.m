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
toSave = 'yes';
%% SET RUN PARAMS
params.alignEvent          = 'goCue';
params.lowFR               = 1; % remove clusters firing less than this val
params.dt = 0.05;
params.moveThresh          = 0.15;

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
%     jawAngle = getJawAngle(taxis, obj, met);         % Find tongue angle

    kin = struct();                                    % Make an empty struct for kinematic information
    conditions = {1,2};

    numTrials = 0;                                     % Get total number of trials that you are working with
    for i = 1:numel(conditions)
        cond = conditions{i};
        nums = length(met.trialid{cond});
        numTrials = numTrials + nums;
    end


    % Find kinematic data for all DeepLabCut features
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

    % Find interpolated ME for each condition
    edges = 0:0.005:5.5;
    [met,mov,me] = assignEarlyTrials(obj,met,params);
    [MEinterp,~] = findInterpME(edges,conditions, met,mov,me);
    kin.MotionEnergy = [MEinterp{1} MEinterp{2}];

    modeparams.tix = 1:1100;       % time points to use when finding mode
    modeparams.fcut = 50;          % smoothing cutoff frequency
    modeparams.cond = 1:2;         % which conditions to use to find mode
    modeparams.method = 'xcorr';   % 'xcorr' or 'regress' (basically the same)
    modeparams.fa = false;         % if true, reduces neural dimensions to 10 with factor analysis

    % get modes based on single trial full neural data (or latents) and kinematic features
    kinfns = fieldnames(kin);
    for i = 1:numel(kinfns)         % For all features...
        Y = kin.(kinfns{i});            % feature data to use to calculate mode
        [mode.(kinfns{i}), dat.(kinfns{i})] = findMode(obj, Y, modeparams);     % Find the mode for the current feature
        proj.(kinfns{i}) = getProjection(dat.(kinfns{i}), mode.(kinfns{i}));    % Project all single trials onto this mode
    end

    % orthogonalize modes to each other
    kinmodes = zeros(numel(mode.(kinfns{1})),numel(kinfns));
    for i = 1:numel(kinfns)
        kinmodes(:,i) = mode.(kinfns{i});
    end
    % TODO: kinmodes has rank of 1, so orthModes has 1 column vector in it
    orthModes = gschmidt(kinmodes);

    for i = 1:numel(kinfns)
        orthmode.(kinfns{i}) = orthModes(:,i);
    end


    % project data onto orthogonalized modes
    for i = 1:numel(kinfns)
        orthproj.(kinfns{i}) = getProjection(dat.(kinfns{i}), orthmode.(kinfns{i}));
    end

    % mode is a struct with unorthgonalized modes
    % orthmode is a struct with orthogonalized modes
    % proj is a struct with single trials projected onto modes
    % orthproj is a struct with single trials projected onto orthmodes

    
    % plot everything
    for i = 1:numel(kinfns)
        f(i) = subplot(3,2,i);
        imagesc(taxis,1:numTrials,orthproj.(kinfns{i})')
        colorbar(f(i))
        xlabel('Time since go-cue (s)')
        title(kinfns{i});
    end

    sesstitle = strcat(anm,date,' ;  ','Probe ',probenum,'Kinematic Modes');  % Name/title for session
    sgtitle(sesstitle,'FontSize',16)

    % next steps:
    % Qpotent = orthModes;
    % Qnull is then the PCs of the residual activity after removing Qpotent
    % from single trials (or PSTHs)

end

