% Script for building linear decoder of trial type based on jaw velocity
% and choide coding dimension
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
%% Load all of the data
sesh = 15;
obj = objs{sesh};     % 11th data object = JEB7, 04-29 (Classic sesh)
met = meta(sesh);

anm = obj.pth.anm;                  % Animal name
date = obj.pth.dt;                  % Session date
probenum = string(met.probe);       % Which probe was used
%% Find jaw vel
% Find the jaw velocity at all time points in the session for trials of
% specific conditions
conditions = {1,2};
if strcmp(params.jawMeasure,'sideJaw')
    jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'vel');
elseif strcmp(params.jawMeasure,'Trident')
    jaw_by_cond = findTridentVelocity(taxis, obj,conditions,met);
end
%% Find choice mode
rez.time = objs{1}.time;
rez.condition = objs{1}.condition;
rez.alignEvent = params.alignEvent;

% Find CDchoice (coding dimension during delay period)
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
epoch = 'latedelay';
choice_mode = choiceMode(obj,met,cond,epoch,rez.alignEvent,'no');

% Project single trials onto choice mode
cd = choice_mode;
latent = getTrialLatents(obj,cd,conditions,met);
lat_choice = [];
jaw = [];
for c = 1:numel(conditions)
    lat_choice = [lat_choice,latent{c}];
    jaw = [jaw,jaw_by_cond{c}];
end
%% Find average jaw vel and choice mode for each trial
% Define time intervals: Time frame for late delay period(from -0.4 before go-cue to -0.1)
late_start = find(taxis>=-0.4, 1, 'first');
late_stop = find(taxis<=-0.05, 1, 'last');
lateDelay = late_start:late_stop;

% Get jaw velocity and activity mode averages for late delay
timeInt = lateDelay;
jawVel_late = getAverages(timeInt,jaw);
Choice_late = getAverages(timeInt,lat_choice);
%%  Format data properly for the SVM function
X1 = jawVel_late;        % Avg jaw velocity during late delay on each trial
X2 = Choice_late;        % Avg choice mode value during late delay on each trial
Y = NaN(1,length(X1));   % Classification of each trial as 'R' or 'L'; 1 = R and 2 = L

cnt = 0;
for c = 1:length(conditions)
    ntrix = length(met.trialid{c});
    if cnt == 0
        Y(1:ntrix) = 1;     % For the first condition, classify each trial as 1 ('Right trials')
    else
        Y(cnt+1:end) = 2;   % For the second condition, classify each trial as 2 ('Left trials')
    end
    cnt = ntrix;
end

nanix = find(isnan(X1));                    % Find indices where the jaw vel is a NaN
X1 = X1(~isnan(X1));                        % Get rid of the NaN values
X2(nanix) = [];                             % Indices that were a NaN for jaw vel, get rid of those indices in the choice mode as well
Y(nanix) = [];                              % Get rid of trials that had NaN values for jaw or choice
%% Define train and test set for linear decoder
trainpct = 0.5;         % Fraction of trials that you want to be training trials for the SVM
testpct = 1-trainpct;   % Fraction of trials that you want to be testing trials for the SVM
totTrials = length(Y);

numTrain = floor(totTrials*trainpct);               % Number of trials for the training set
trainix = sort(randsample(totTrials,numTrain));     % Randomly sample trials for the training set

X1_train = (X1(trainix))';          % Jaw training set
X2_train = (X2(trainix))';          % Choice mode training set
Y_train = (Y(trainix))';            % Trial type for training set 

test = ones(1,totTrials);       
test(trainix) = 0;
testix = find(test)';               % Use the remaining trials as the test set

X1_test = (X1(testix))';
X2_test = (X2(testix))';
Y_test = (Y(testix))';

% Fit model on training data
model_jaw = fitclinear(X1_train,Y_train);
model_choice = fitclinear(X2_train,Y_train);

% Test model on testing data
[label.jaw,score.jaw] = predict(model_jaw,X1_test);
[label.choice,score.choice] = predict(model_choice,X2_test);

% Compare model predictions to ground truth
numCorrect.jaw = label.jaw==Y_test;
numCorrect.choice = label.choice==Y_test;
pctCorrect.jaw = sum(numCorrect.jaw)/length(Y_test);
pctCorrect.choice = sum(numCorrect.choice)/length(Y_test);

%%%%% How to get a 'score' for each trial? So that I can set a threshold
%%%%% and generate an ROC curve for a single session?

disp(pctCorrect)