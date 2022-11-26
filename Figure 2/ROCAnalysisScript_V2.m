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
AUC.jaw = NaN(1,length(meta));          % (1 x num sessions)
AUC.choice = NaN(1,length(meta));       % (1 x num sessions)
R = NaN(1,length(meta));                % Store R^2 values for jaw vel vs choice for each session
for gg = 1:length(meta)
    figure(gg);
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
    lat_choice = [];
    jaw = [];
    for c = 1:numel(conditions)
        lat_choice = [lat_choice,latent{c}];
        jaw = [jaw,jaw_by_cond{c}];
    end
    %%%% Find average jaw vel and choice mode for each trial %%%%
    % Define time intervals: Time frame for late delay period(from -0.4 before go-cue to -0.1)
    late_start = find(taxis>=-0.4, 1, 'first');
    late_stop = find(taxis<=-0.05, 1, 'last');
    lateDelay = late_start:late_stop;

    % Get jaw velocity and activity mode averages for late delay
    timeInt = lateDelay;
    jawVel_late = getAverages(timeInt,jaw);
    Choice_late = getAverages(timeInt,lat_choice);

    %%%%  Format data properly for the function %%%%
    jv = jawVel_late;        % Avg jaw velocity during late delay on each trial
    ch = Choice_late;        % Avg choice mode value during late delay on each trial
    Y = NaN(1,length(jv));   % Classification of each trial as 'R' or 'L'; 0 = R and 1 = L

    cnt = 0;
    for c = 1:length(conditions)
        ntrix = length(met.trialid{c});
        if cnt == 0
            Y(1:ntrix) = 0;     % For the first condition, classify each trial as 0 ('Right trials')
        else
            Y(cnt+1:end) = 1;   % For the second condition, classify each trial as 1 ('Left trials')
        end
        cnt = ntrix;
    end

    nanix = find(isnan(jv));                    % Find indices where the jaw vel is a NaN
    jv = jv(~isnan(jv));                        % Get rid of the NaN values
    ch(nanix) = [];                             % Indices that were a NaN for jaw vel, get rid of those indices in the choice mode as well
    Y(nanix) = [];                              % Get rid of trials that had NaN values for jaw or choice

    R(gg) = corr2(jv,ch);        % Get R^2 value for current session (how correlated jaw vel is to choice mode on a trial by trial basis)

    %%%% DO ROC analysis for jaw velocity %%%%
    for oo = 1:2                                            % For jaw velocity and choice mode...
        if oo==1
            metric = jv;
            col = 'cyan';
        elseif oo==2
            metric = ch;
            figtitle = 'ROC analysis for choice mode';
            col = 'magenta';
        end
        nsteps = 10000;                                      % How many steps you want to do the ROC analysis for
        rang = linspace(min(metric),max(metric),nsteps-2);   % Generate steps from the min value of jaw vel or choice mode
        dr = rang(2)-rang(1);
        steps = NaN(1,nsteps);
        steps(2:end-1) = rang; steps(1) = steps(2)-dr; steps(end) = steps(end-1)+dr;   % Add one more step at the beginning and end of the range


        trueposrate = NaN(1,nsteps);                   % Store true positive and false positive rates for each threshold step
        falseposrate = NaN(1,nsteps);
        for i = 1:nsteps
            thresh = steps(i);                         % Choose the threshold
            posix = find(Y==1);
            posavg = mean(metric(posix));
            negix = find(Y==0);
            negavg = mean(metric(negix));
            if posavg>negavg                               % If the maximum value of jaw vel or choice is a positive class (left trial), want to classify any values that are greater than thresh as positive
                prediction = metric>thresh;                % Greater than or equal to threshold = 1; Less than or equal to threshold = 0;
            else                                           % If the minimum value of jaw vel or choice is a positive class (left trial), want to classify any values that are less than thresh as positive
                prediction = metric<thresh;                % Less than or equal to threshold = 1; Greater than or equal to threshold = 0;
            end
            corrpos = sum(prediction==1 & Y==1);       % Positives correctly classified (when they are predicted as 1 and they really are 1)
            incorrneg = sum(prediction==1 & Y==0);     % Negatives incorrectly classified (when they are predicted as 1 but they are really negative)
            totalpos = sum(Y);                         % Total num trials that are truly positive
            totalneg = length(Y)-totalpos;             % Total num trials that are truly negative
            trueposrate(i) = corrpos/totalpos;         % True pos rate = positives correctly classified/total true positives
            falseposrate(i) = incorrneg/totalneg;      % False pos rate = negatives incorrectly classified/total true negatives
        end

        if issorted(falseposrate,'ascend')                 % Ordering of the x-axis values matters for AUC magnitude
            auc = trapz(falseposrate,trueposrate);         % Find area under the curve for the ROC plot
        else
            auc = trapz(flip(falseposrate),flip(trueposrate));
        end

        if oo ==1
            AUC.jaw(gg) = auc;          % Save AUC to correct metric and session
            strauc = num2str(auc);
            if strcmp(params.jawMeasure,'sideJaw')
                figtitle1 = strcat('Jaw velocity; AUC =',' ',strauc);
            elseif strcmp(params.jawMeasure,'MotionEnergy')
                figtitle1 = strcat('Motion Energy; AUC =',' ',strauc);
            end
        elseif oo==2
            AUC.choice(gg) = auc;
            strauc = num2str(auc);
            figtitle2 = strcat('CDlate; AUC =',' ',strauc);
        end

        %%% PLOT ROC CURVES FOR GIVEN SESSION %%%
        figure(gg);
        x = linspace(0,1);
        y = x;
        plot(falseposrate,trueposrate,'LineWidth',2,'Color',col); hold on;
        xlabel('False positive rate','FontSize',13)
        ylabel('True positive rate','FontSize',13)
        blah = strcat('Example ROC Analysis for',anm,date,'Probe',probenum);
        sgtitle(blah)
    end
    plot(x,y,'LineStyle','--','Color','black');
    legend(figtitle1,figtitle2,'FontSize',12,'Location','best')
end
%% Visualize AUC stats across sessions
figure();
scatter(AUC.jaw,AUC.choice,'black','filled'); hold on;
scatter(AUC.jaw(3),AUC.choice(3),'green','filled')
if strcmp(params.jawMeasure,'sideJaw')
    xlabel('AUC for jaw velocity','FontSize',13);
elseif strcmp(params.jawMeasure,'MotionEnergy')
    xlabel('AUC for Motion Energy','FontSize',13);
end
ylabel('AUC for CDlate','FontSize',13);
xlim([0.45 1])
xline(0.5,'LineStyle','--')
yline(0.5,'LineStyle','--')
ylim([0.45 1])
title('Trial type decoding across all sessions; avg from late delay','FontSize',14)

goodjaw = sum(AUC.jaw>0.7)/length(AUC.jaw);
goodchoice = sum(AUC.choice>0.7)/length(AUC.choice);
disp(goodjaw)
disp(goodchoice)
%% Histogram of R^2 values across sessions
figure();
nbins = 12;
R = abs(R);
histogram(R,nbins,'FaceColor','black','FaceAlpha',0.5)
xlabel('R^2 absolute value','FontSize',13)
ylabel('# sessions','FontSize',13)
if strcmp(params.jawMeasure,'sideJaw')
    title('Correlation between jaw velocity and CDlate','FontSize',14)
elseif strcmp(params.jawMeasure,'MotionEnergy')
    title('Correlation between motion energy and CDlate','FontSize',14)
end
