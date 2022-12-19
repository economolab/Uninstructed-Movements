% Find AVG CDs across all sessions
% Find AVG movement across each trial type
% Example session of each

clear,clc,close all

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'
params.jawMeasure          = 'MotionEnergy'; % sideJaw or Trident
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

%% SET METADATA

meta = [];
meta = loadJEB6_ALMVideo(meta);
meta = loadJEB7_ALMVideo(meta);
meta = loadEKH1_ALMVideo(meta);
meta = loadEKH3_ALMVideo(meta);
meta = loadJGR2_ALMVideo(meta);
meta = loadJGR3_ALMVideo(meta);
meta = loadJEB15_ALMVideo(meta);
meta = loadJEB14_ALMVideo(meta);
% meta = loadJEB14_M1TJVideo(meta);
% meta = loadEKH3_M1TJVideo(meta);

params.probe = [meta.probe]; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

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
inclJEB15 = 'yes';
if strcmp(inclJEB15,'yes')
    anm = 'JEB15';
    dates = {'2022-07-26','2022-07-27','2022-07-28'};       % Dates that you want to combine probe info from
    [meta, objs, params] = combineSessionProbes(meta,objs,params,anm,dates);
end
%% Find average ME or jaw velocity for each trial type
for oo = 1:length(meta)
    obj = objs{oo};
    met = meta(oo);
    % Find the jaw velocity at all time points in the session for trials of
    % specific conditions
    conditions = {1,2};
    taxis = obj.time;
    if strcmp(params.jawMeasure,'sideJaw')
        jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'vel',params);
    elseif strcmp(params.jawMeasure,'MotionEnergy')
        params.moveThresh          = 0.15;      % What percentage of the delay period you want to use for identifying early move trials
        [met,mov,me] = assignEarlyTrials(obj,met,params);
        [jaw_by_cond,~] = findInterpME(taxis,conditions, met,mov,me,params,obj);
    end
    moveselect = NaN(length(taxis),2);
    moveselect(:,1) = mean(jaw_by_cond{1},2); moveselect(:,2) = mean(jaw_by_cond{2},2);
    rez(oo).moveselect = moveselect;
end
%% Find all CDs for each session
for sesh = 1:length (meta)
    obj = objs{sesh};
    met = meta(sesh);

    anm = obj.pth.anm;                  % Animal name
    date = obj.pth.dt;                  % Session date
    probenum = string(met.probe);       % Which probe was used

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
end
%% Avg trial-type selective movements and CDlate across sessions (and get standard deviation)
tempmove.R = [];
tempmove.L = [];
tempCD.late.R = []; tempCD.early.R = []; tempCD.go.R = [];
tempCD.late.L = []; tempCD.early.L = []; tempCD.go.L = [];
for oo = 1:length(meta)
    tempmove.R = [tempmove.R, rez(oo).moveselect(:,1)];  tempmove.L = [tempmove.L, rez(oo).moveselect(:,2)];
    tempCD.late.R = [tempCD.late.R, rez(oo).cdLate_latent(:,1)];  tempCD.late.L = [tempCD.late.L, rez(oo).cdLate_latent(:,2)];
    tempCD.early.R = [tempCD.early.R, rez(oo).cdEarly_latent(:,1)];  tempCD.early.L = [tempCD.early.L, rez(oo).cdEarly_latent(:,2)];
    tempCD.go.R = [tempCD.go.R, rez(oo).cdGo_latent(:,1)];  tempCD.go.L = [tempCD.go.L, rez(oo).cdGo_latent(:,2)];
end
movement.avg = NaN(length(taxis),2);
movement.avg(:,1) = mean(tempmove.R,2,'omitnan'); movement.avg(:,2) = mean(tempmove.L,2,'omitnan');

cdlate.avg = NaN(length(taxis),2);  cdlate.std = NaN(length(taxis),2);
cdlate.avg(:,1) = mean(tempCD.late.R,2,'omitnan'); cdlate.avg(:,2) = mean(tempCD.late.L,2,'omitnan');
cdlate.std(:,1) = std(tempCD.late.R,0,2,'omitnan'); cdlate.std(:,2) = std(tempCD.late.L,0,2,'omitnan');

cdearly.avg = NaN(length(taxis),2);  cdearly.std = NaN(length(taxis),2);
cdearly.avg(:,1) = mean(tempCD.early.R,2,'omitnan'); cdearly.avg(:,2) = mean(tempCD.early.L,2,'omitnan');
cdearly.std(:,1) = std(tempCD.early.R,0,2,'omitnan'); cdearly.std(:,2) = std(tempCD.early.L,0,2,'omitnan');

cdgo.avg = NaN(length(taxis),2); cdgo.std = NaN(length(taxis),2);
cdgo.avg(:,1) = mean(tempCD.go.R,2,'omitnan'); cdgo.avg(:,2) = mean(tempCD.go.L,2,'omitnan');
cdgo.std(:,1) = std(tempCD.go.R,0,2,'omitnan'); cdgo.std(:,2) = std(tempCD.go.L,0,2,'omitnan');
%% Get confidence intervals for all CDs
nSessions = length(meta);
upperci.early = NaN(length(taxis)-1,2); lowerci.early = NaN(length(taxis)-1,2);
upperci.late = NaN(length(taxis)-1,2); lowerci.late = NaN(length(taxis)-1,2);
upperci.go = NaN(length(taxis)-1,2); lowerci.go = NaN(length(taxis)-1,2);
for c = 1:numel(conditions)
    upperci.early(:,c) = cdearly.avg(2:end,c)+1.96*(cdearly.std(2:end,c)/nSessions);  % Find the upper 95% confidence interval for each condition
    lowerci.early(:,c) = cdearly.avg(2:end,c)-1.96*(cdearly.std(2:end,c)/nSessions);  % Find lower 95% condifence interval for each condition

    upperci.late(:,c) = cdlate.avg(2:end,c)+1.96*(cdlate.std(2:end,c)/nSessions);  % Find the upper 95% confidence interval for each condition
    lowerci.late(:,c) = cdlate.avg(2:end,c)-1.96*(cdlate.std(2:end,c)/nSessions);  % Find lower 95% condifence interval for each condition

    upperci.go(:,c) = cdgo.avg(2:end,c)+1.96*(cdgo.std(2:end,c)/nSessions);  % Find the upper 95% confidence interval for each condition
    lowerci.go(:,c) = cdgo.avg(2:end,c)-1.96*(cdgo.std(2:end,c)/nSessions);  % Find lower 95% condifence interval for each condition

end
%% Plot avg CDs across all sessions
colors = {[0 0.4470 0.7410],[0.6350 0.0780 0.1840]};
figure();
subplot(3,1,1)
for c = 1:numel(conditions)
    plot(taxis,cdearly.avg(:,c),'LineWidth',2); hold on; 
    patch([taxis(10:end) fliplr(taxis(10:end))],[lowerci.early(9:end,c)' fliplr(upperci.early(9:end,c)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
    ylim([-100 0])
end
title('CDEarly')
legend('Right', 'R CI','Left','L CI')

subplot(3,1,2)
for c = 1:numel(conditions)
    plot(taxis,cdlate.avg(:,c),'LineWidth',2); hold on; 
    patch([taxis(10:end) fliplr(taxis(10:end))],[lowerci.late(9:end,c)' fliplr(upperci.late(9:end,c)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
ylim([-80 20])
end
title('CDLate')
legend('Right', 'R CI','Left','L CI')

subplot(3,1,3)
for c = 1:numel(conditions)
    plot(taxis,cdgo.avg(:,c),'LineWidth',2); hold on; 
    patch([taxis(10:end) fliplr(taxis(10:end))],[lowerci.go(9:end,c)' fliplr(upperci.go(9:end,c)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
ylim([-60 60])
end
title('CDGo')
legend('Right', 'R CI','Left','L CI')
