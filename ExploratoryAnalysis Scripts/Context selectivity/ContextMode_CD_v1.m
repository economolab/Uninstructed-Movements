% Find AVG context mode across all sessions
% Find AVG movement across each context type
% Example session of each

clear,clc,close all

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\ExploratoryAnalysis Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\cd_code_jackie'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\functions'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\utils'));
%% SET RUN PARAMS
params.alignEvent          = 'firstLick'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'
params.jawMeasure          = 'sideJaw'; % sideJaw or Trident
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val
params.moveThresh          = 0.15;

% set conditions to calculate PSTHs for
params.condition(1)     = {'hit&~stim.enable&~autowater&~early'};         % all 2AFC hits, no stim
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};          % all AW hits, no stim
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};        
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         

% set conditions used for finding activity modes
% aw = '2'; % 1-on, 2-off
% stim = '0'; % 0-off
% params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %rhit, aw off
% params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %lhit, aw off
% params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %rmiss, aw off
% params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %lmiss, aw off
% params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};    % hit, aw off
% params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};    % hit, aw off

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
% meta = loadJEB15_ALMVideo(meta);
% meta = loadJEB14_ALMVideo(meta);
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
inclJEB15 = 'no';
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
    moveselect(:,1) = mean(jaw_by_cond{1},2,'omitnan'); moveselect(:,2) = mean(jaw_by_cond{2},2,'omitnan');
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

%     % Calculate all CDs
%     cond = [3,4];
%     [Early, Late, Go] = findAllAFCCDs(rez,sesh,ev,params,cond);
%     rez(sesh).cdEarly_mode = Early;
%     rez(sesh).cdLate_mode = Late;  rez(sesh).cdGo_mode = Go;

    % Calculate Context mode
    cond = [1,2];
    e1 = mode(ev.sample) - 0.3 - median(ev.(params.alignEvent));
    e2 = mode(ev.sample) - 0.05 - median(ev.(params.alignEvent));

    times.context = rez(sesh).time>e1 & rez(sesh).time<e2;
    cdContext_mode = calcAFCCD(rez(sesh),times.context,cond);
    rez(sesh).cdContext_mode = cdContext_mode;


    % orthogonalize
    [fns,~] = patternMatchCellArray(fieldnames(rez(sesh)),{'mode'},'all');
    modes = zeros(numel(rez(sesh).cdContext_mode),numel(fns));
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
%% Plot Context Mode and Avg trial-type movements
colors = {[0.25 0.25 0.25],'magenta'};
for sesh = 1:length(meta)
    met = meta(sesh);
    obj = objs{sesh};
    plotit = strcat(met.anm,{' '},met.date);
    
    figure();
    subplot(2,1,1)
    for c = 1:numel(conditions)
        plot(taxis,rez(sesh).cdContext_latent(:,c),'Color',colors{c},'LineWidth',2)
        hold on;
    end
    legend('2AFC','AW')
    title('Context Mode')

    subplot(2,1,2)
    for c = 1:numel(conditions)
        plot(taxis,rez(sesh).moveselect(:,c),'Color',colors{c},'LineWidth',2)
        hold on;
    end
    legend('2AFC','AW')
    if strcmp(params.jawMeasure,'sideJaw')
        title('Avg jaw vel')
    else
        title('Avg Motion Energy')
    end

    sgtitle(plotit)
end

%% Avg trial-type selective movements and CDlate across sessions (and get standard deviation)
tempmove.AFC = [];
tempmove.AW = [];
tempCD.AFC = [];
tempCD.AW = [];
for oo = 1:length(meta)
    tempmove.AFC = [tempmove.AFC, rez(oo).moveselect(:,1)];  tempmove.AW = [tempmove.AW, rez(oo).moveselect(:,2)];
    tempCD.AFC = [tempCD.AFC, rez(oo).cdContext_latent(:,1)];  tempCD.AW = [tempCD.AW, rez(oo).cdContext_latent(:,2)];
end
movement.avg = NaN(length(taxis),2);  movement.std = NaN(length(taxis),2);
movement.avg(:,1) = mean(tempmove.AFC,2,'omitnan'); movement.avg(:,2) = mean(tempmove.AW,2,'omitnan');
movement.std(:,1) = std(tempmove.AFC,0,2,'omitnan'); movement.std(:,2) = std(tempmove.AW,0,2,'omitnan');

cdcontext.avg = NaN(length(taxis),2);  cdcontext.std = NaN(length(taxis),2);
cdcontext.avg(:,1) = mean(tempCD.AFC,2,'omitnan'); cdcontext.avg(:,2) = mean(tempCD.AW,2,'omitnan');
cdcontext.std(:,1) = std(tempCD.AFC,0,2,'omitnan'); cdcontext.std(:,2) = std(tempCD.AW,0,2,'omitnan');

%% Get confidence intervals for all CDs
nSessions = length(meta);
upperci.context = NaN(length(taxis)-1,2); lowerci.context = NaN(length(taxis)-1,2);
upperci.move = NaN(length(taxis)-1,2); lowerci.move = NaN(length(taxis)-1,2);

for c = 1:numel(conditions)
    upperci.context(:,c) = cdcontext.avg(2:end,c)+1.96*(cdcontext.std(2:end,c)/nSessions);  % Find the upper 95% confidence interval for each condition
    lowerci.context(:,c) = cdcontext.avg(2:end,c)-1.96*(cdcontext.std(2:end,c)/nSessions);  % Find lower 95% condifence interval for each condition

    upperci.move(:,c) = movement.avg(2:end,c)+1.96*(movement.std(2:end,c)/nSessions);  % Find the upper 95% confidence interval for each condition
    lowerci.move(:,c) = movement.avg(2:end,c)-1.96*(movement.std(2:end,c)/nSessions);  %
end
%% Plot avg CDs across all sessions
colors = {[0.25 0.25 0.25],'magenta'};
figure();
subplot(2,1,1)
for c = 1:numel(conditions)
    plot(taxis,cdcontext.avg(:,c),'LineWidth',2); hold on;
    patch([taxis(10:end) fliplr(taxis(10:end))],[lowerci.context(9:end,c)' fliplr(upperci.context(9:end,c)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
end
title('CDContext')
legend('2AFC', '2AFC CI','AW','AW CI')

subplot(2,1,2)
for c = 1:numel(conditions)
    plot(taxis,movement.avg(:,c),'LineWidth',2); hold on;
    patch([taxis(10:end) fliplr(taxis(10:end))],[lowerci.move(9:end,c)' fliplr(upperci.move(9:end,c)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
end
title('Movements')
legend('2AFC', '2AFC CI','AW','AW CI')

%%

% Munib's latest version for finding CDearly, late, and go
function [cdEarly_mode, cdLate_mode, cdGo_mode] = findAllAFCCDs(rez,sessix,ev,params,cond)

% cd early mode
e1 = mode(ev.sample) + 0.4 - mode(ev.(params.alignEvent));
e2 = mode(ev.sample) + 0.8 - mode(ev.(params.alignEvent));

times.early = rez(sessix).time>e1 & rez(sessix).time<e2;
cdEarly_mode = calcAFCCD(rez(sessix),times.early,cond);


% cd late mode
e1 = mode(ev.(params.alignEvent)) - 0.6 - mode(ev.(params.alignEvent));
e2 = mode(ev.(params.alignEvent)) - 0.15 - mode(ev.(params.alignEvent));

times.late = rez(sessix).time>e1 & rez(sessix).time<e2;
cdLate_mode = calcAFCCD(rez(sessix),times.late,cond);


% cd go mode
e1 = mode(ev.(params.alignEvent)) + 0.02 - mode(ev.(params.alignEvent));
e2 = mode(ev.(params.alignEvent)) + 0.42 - mode(ev.(params.alignEvent));

times.go = rez(sessix).time>e1 & rez(sessix).time<e2;
cdGo_mode = calcAFCCD(rez(sessix),times.go,cond);
end

function cd = calcAFCCD(rez,times,cond)
tempdat = rez.psth(:,:,cond);
mu = squeeze(mean(tempdat(times,:,:),1));
sd = squeeze(std(tempdat(times,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
end

