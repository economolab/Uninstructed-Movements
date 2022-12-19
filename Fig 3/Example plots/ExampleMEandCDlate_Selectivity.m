clear,clc,close all

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Figure2'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\functions'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Utils'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\cd_code_jackie'));
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
% meta = loadJEB6_ALMVideo(meta);
meta = loadJEB7_ALMVideo(meta);
% meta = loadEKH1_ALMVideo(meta);
% meta = loadEKH3_ALMVideo(meta);
% meta = loadJGR2_ALMVideo(meta);
% meta = loadJGR3_ALMVideo(meta);
% meta = loadJEB15_ALMVideo(meta);
% meta = loadJEB14_ALMVideo(meta);

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
% anm = 'JEB15';
% dates = {'2022-07-26','2022-07-27','2022-07-28'};       % Dates that you want to combine probe info from
% [meta, objs, params] = combineSessionProbes(meta,objs,params,anm,dates);
%% Find Motion Energy for all trials in example session
sesh = 2;           % Maybe 1,4
obj = objs{sesh};     
met = meta(sesh);

anm = obj.pth.anm;                  % Animal name
date = obj.pth.dt;                  % Session date
probenum = string(met.probe);       % Which probe was used

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
l1 = size(jaw_by_cond{1},2);      % Number of trials in the first condition

%% Find CDlate for example session
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
%% Get single trial projections onto CDlate

obj = objs{sesh};     % 11th data object = JEB7, 04-29 (Classic sesh)
met = meta(sesh);

anm = obj.pth.anm;                  % Animal name
date = obj.pth.dt;                  % Session date
probenum = string(met.probe);       % Which probe was used

rez(sesh).time = objs{1}.time;
rez(sesh).alignEvent = params.alignEvent;

% Project single trials onto choice mode
cd = rez(sesh).cdLate_mode;
latent = getTrialLatents(obj,cd,conditions,met);        % (1 x num conditions); same dimensions as 'jaw_by_cond'
%% Get averages and standard deviation for ME and CDlatents
for c = 1:numel(cond)
    avg.cdLate_latent(:,c) = mySmooth(mean(latent{c},2),51);
    avg.MEinterp(:,c) = mySmooth(mean(jaw_by_cond{c},2),51);

    stdev.cdLate_latent(:,c) = std(latent{c},0,2);
    stdev.MEinterp(:,c) = std(jaw_by_cond{c},0,2);
end
%%
alpha = 0.1; % transparency of shaded confidence interval / std error

% Plot
figure();
colors = {[0 0.4470 0.7410],[0.6350 0.0780 0.1840]};
for c = 1:numel(cond)
    temp = cond(c);
    nTrials = size(jaw_by_cond{temp},2);            % Needed to calculate standard error of mean for error bars
    col = colors{c};
    subplot(2,1,2)
    plot(taxis,avg.cdLate_latent(:,temp),'Color',col,'LineWidth',3)
    hold on;
    ax = gca;
    shadedErrorBar(taxis, avg.cdLate_latent(:,temp), stdev.cdLate_latent(:,temp) ./ sqrt(nTrials) ,{'Color',col,'LineWidth',2}, alpha, ax)
    set(gca, 'YDir','reverse')
    ylabel('a.u.')
    xlim([-2.3 0])
    xlabel('Time since go-cue')
    
    subplot(2,1,1)
    plot(taxis,avg.MEinterp(:,temp),'Color',col,'LineWidth',2.5); hold on;
    ax = gca;
    shadedErrorBar(taxis, avg.MEinterp(:,temp), stdev.MEinterp(:,temp) ./ sqrt(nTrials) ,{'Color',col,'LineWidth',2}, alpha, ax)
    ylabel('M.E.')
    xlim([-2.3 0])
    xlabel('Time since go-cue')
end


legend('Right','Left','Location','best')