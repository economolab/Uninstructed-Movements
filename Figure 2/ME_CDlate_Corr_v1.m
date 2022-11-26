% Script for doing ROC analysis--ability to decode trial type from jaw velocity and choice
% mode
clear; clc; close all;
%%
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\ActivityModes'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Utils'));
%% SET RUN PARAMS

params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'
params.jawMeasure          = 'MotionEnergy'; % sideJaw or MotionEnergy
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val
params.moveThresh          = 0.15;

% set conditions to calculate PSTHs for
params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off

% set conditions used for finding activity modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %rhit, aw off
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %lhit, aw off
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %rmiss, aw off
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %lmiss, aw off
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};    % hit, aw off

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 31;

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'excellent','great','good','multi','fair','poor'};

%% SET METADATA FROM ALL RELEVANT SESSIONS/ANIMALS
meta = [];
meta = loadJEB6_ALMVideo(meta);
meta = loadJEB7_ALMVideo(meta);
meta = loadEKH1_ALMVideo(meta);
meta = loadEKH3_ALMVideo(meta);
meta = loadJGR2_ALMVideo(meta);
meta = loadJGR3_ALMVideo(meta);
meta = loadJEB14_ALMVideo(meta);
meta = loadJEB15_ALMVideo(meta);

params.probe = [meta.probe];

taxis = meta(end).tmin:meta(end).dt:meta(end).tmax;   % get time-axis with 0 as time of event you aligned to
taxis = taxis(1:end-1);
%% LOAD AND PROCESS DATA

objs = loadObjs(meta);

for metaix = 1:numel(meta)
    obj = objs{metaix};
    disp('______________________________________________________')
    disp(['Processing data for session ' [meta(metaix).anm '_' meta(metaix).date]])             % Display progress
    disp(' ')
    [sessparams{metaix},sessobj{metaix}] = processData(obj,params,params.probe(metaix));        % Function for processing all of the data objs
end

% clean up sessparams and sessobj
for metaix = 1:numel(meta)
    params.trialid{metaix} = sessparams{metaix}.trialid;
    params.cluid{metaix} = sessparams{metaix}.cluid{params.probe(metaix)};

    objs{metaix} = sessobj{metaix};
    objs{metaix}.psth = objs{metaix}.psth{params.probe(metaix)};
    objs{metaix}.trialdat = objs{metaix}.trialdat{params.probe(metaix)};
    objs{metaix}.presampleFR = objs{metaix}.presampleFR{params.probe(metaix)};
    objs{metaix}.presampleSigma = objs{metaix}.presampleSigma{params.probe(metaix)};
end

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')
%%
% Remove unwanted sessions
[meta,objs,params] = useInclusionCriteria(objs,params,meta);
%% Adjust for older functions
for i = 1:numel(objs)
    meta(i).trialid = params.trialid{i};
    conditions = {1,2};

    for c = 1:numel(conditions)
        trix = meta(i).trialid{c};
        objs{i}.trialpsth_cond{c} = objs{i}.trialdat(:,:,trix);
    end
end
%% Combine cells from both probes in a session
anm = 'JEB15';
dates = {'2022-07-26','2022-07-27','2022-07-28'};       % Dates that you want to combine probe info from
[meta, objs, params] = combineSessionProbes(meta,objs,params,anm,dates);
%% Load all of the data
R.left = NaN(1,length(meta));                % Store R^2 values for jaw vel vs choice for each session
R.right = NaN(1,length(meta));
for gg = 1:length(meta)
    sesh = gg;
    obj = objs{sesh};     % 11th data object = JEB7, 04-29 (Classic sesh)
    met = meta(sesh);

    anm = obj.pth.anm;                  % Animal name
    date = obj.pth.dt;                  % Session date
    probenum = string(met.probe);       % Which probe was used

    %%%% FIND JAW VEL %%%%
    % Find the jaw velocity at all time points in the session for trials of
    % specific conditions
    conditions = {1,2};
    if strcmp(params.jawMeasure,'sideJaw')
        jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'vel',params);
    elseif strcmp(params.jawMeasure,'MotionEnergy')
        [met,mov,me] = assignEarlyTrials(obj,met,params);
        jaw_by_cond = findInterpME(taxis,conditions, met,mov,me,params,obj);
    end

    %%%% FIND CHOICE MODE %%%%
    % Find CDchoice (coding dimension during delay period)
    ev.sample = objs{sesh}.bp.ev.sample;
    ev.delay = objs{sesh}.bp.ev.delay;
    ev.goCue = objs{sesh}.bp.ev.goCue;
    ev.(params.alignEvent) = objs{sesh}.bp.ev.(params.alignEvent);

    rez(sesh).time = objs{sesh}.time;
    rez(sesh).psth = objs{sesh}.psth;
    rez(sesh).psth = standardizePSTH(objs{sesh});
    rez(sesh).condition = params.condition;
    rez(sesh).alignEvent = params.alignEvent;
    rez(sesh).ev = ev;

    % Calculate all CDs
    [Early, Late, Go] = findAllCDs(rez,sesh,ev,params);
    rez(sesh).cdEarly_mode = Early;
    rez(sesh).cdLate_mode = Late;  rez(sesh).cdGo_mode = Go;

    % orthogonalize
    [fns,~] = patternMatchCellArray(fieldnames(rez(sesh)),{'mode'},'all');
    modes = zeros(numel(rez(sesh).cdLate_mode),numel(fns));
    for i = 1:numel(fns)
        modes(:,i) = rez(sesh).(fns{i});
    end

    orthModes = gschmidt(modes);

    for i = 1:numel(fns)
        rez(sesh).(fns{i}) = orthModes(:,i);
    end

    cond = [1 2];
    for i = 1:numel(fns)
        tempmode = rez(sesh).(fns{i});
        for j = 1:numel(cond)
            c = cond(j);
            tempdat = rez(sesh).psth(:,:,c)*rez(sesh).(fns{i});
            normfactor = 1;
            rez(sesh).([fns{i}(1:end-5) '_latent'])(:,j) = tempdat ./ normfactor;
        end
    end
    % Project single trials onto choice mode
    cd = rez(sesh).cdLate_mode;
    latent = getTrialLatents(obj,cd,conditions,met);
    lat_choice_cond = latent;

    %%%% Find average jaw vel and choice mode for each trial %%%%
    % Define time intervals: Time frame for late delay period(from -0.4 before go-cue to -0.1)
    late_start = find(taxis>=-0.4, 1, 'first');
    late_stop = find(taxis<=-0.05, 1, 'last');
    lateDelay = late_start:late_stop;

    % Get jaw velocity and activity mode averages for late delay
    timeInt = lateDelay;
    for c = 1:numel(conditions)
        tempjaw = jaw_by_cond{c};
        tempcd = lat_choice_cond{c};
        jv{c} = getAverages(timeInt,tempjaw);
        ch{c} = getAverages(timeInt,tempcd);
    end

    for c = 1:numel(conditions)
        jv_T = jv{c}; ch_T = ch{c};
        nanix = find(isnan(jv_T));                    % Find indices where the jaw vel is a NaN
        jv_T = jv_T(~isnan(jv_T));                    % Get rid of the NaN values
        ch_T(nanix) = [];                             % Indices that were a NaN for jaw vel, get rid of those indices in the choice mode as well
        
        if c==1
            R.right(gg) = corr2(jv_T,ch_T);        % Get R^2 value for current session (how correlated jaw vel is to choice mode on a trial by trial basis)
        else
            R.left(gg) = corr2(jv_T,ch_T);
        end
    end
end
%% Histogram of R^2 values across sessions
figure();
nbins = 12;
for c = 1:numel(conditions)
    subplot(2,1,c)
    if c==1
        tempR = abs(R.right);
        col = 'blue';
        plotit = 'Right trials';
    else
        tempR = abs(R.left);
        col = 'red';
        plotit = 'Left trials';
    end
    histogram(tempR,nbins,'FaceColor',col,'FaceAlpha',0.5)
    xlabel('R^2 absolute value','FontSize',13)
    ylabel('# sessions','FontSize',13)
    title(plotit)
end

if strcmp(params.jawMeasure,'sideJaw')
    sgtitle('Correlation between jaw velocity and CDlate','FontSize',14)
elseif strcmp(params.jawMeasure,'MotionEnergy')
    sgtitle('Correlation between motion energy and CDlate','FontSize',14)
end
