% Script for looking at Hazarded Delay Recording + Video
clear,clc,close all

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 31;

params.delay(1) = 0.3000;
params.delay(2) = 0.6000;
params.delay(3) = 1.2000;
params.delay(4) = 1.8000;
params.delay(5) = 2.4000;
params.delay(6) = 3.6000;

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'excellent','great','good','multi','fair','poor'};

%% SET METADATA

meta = [];
meta = loadJEB11_ALMVideo(meta);
meta = loadJEB12_ALMVideo(meta);

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

% remove sessions if:
% 1) less than 40 trials of rhit and lhit each (same as
% 2) atleast 10 clusters of quality that is not noise
use = false(size(objs));
for i = 1:numel(use)
    check(1) = numel(params.trialid{i}{1}) > 40;
    check(2) = numel(params.trialid{i}{2}) > 40;
    check(3) = numel(params.cluid{i}) >= 10;
    if all(check)
        use(i) = true;
    end
end

meta = meta(use);
objs = objs(use);
params.probe = params.probe(use);
params.trialid = params.trialid(use);
params.cluid = params.cluid(use);
%% Adjust for older functions
for i = 1:numel(objs)
    meta(i).trialid = params.trialid{i};
    conditions = {1,2};

    for c = 1:numel(conditions)
        trix = meta(i).trialid{c};
        objs{i}.trialpsth_cond{c} = objs{i}.trialdat(:,:,trix);
    end
end
%% For each session--Trials separated by delay length
plotIndiv = 'no';
for gg = 1:length(meta)
    sesh = gg;
    obj = objs{sesh};
    met = meta(sesh);

    anm = obj.pth.anm;                  % Animal name
    date = obj.pth.dt;                  % Session date
    probenum = string(met.probe);       % Which probe was used

    delaylen = obj.bp.ev.goCue - obj.bp.ev.delay;       % Find the delay length for all trials
    conditions = {1,2};
    met = getDelayTrialID(met,conditions,delaylen);     % Group the trials in each condition based on their delay length
    delPSTH = getPSTHbyDel(params,met,obj);             % Get avg PSTH for each delay length

    % Find the probability of jaw [Jaw] movement at all time points in the session for trials of
    % specific conditions
    taxis = obj.time;
    jaw_by_cond = findJawVelocity(taxis, obj,conditions,met,'prob',params);    % (1 x conditions cell array)
    % Each cell: (time x trials in that condition)

    % Find average jaw velocity for each delay length
    jawvel.left = cell(1,length(params.delay));         % (1 x number of delay lengths)
    jawvel.right = cell(1,length(params.delay));
    for g = 1:length(params.delay)                  % For each delay length...
        gix = find(met.del_trialid{1}==g);              % Get the trial IDs in the first condition that have the current delay length
        tempjaw = nanmean(jaw_by_cond{1}(:,gix),2);     % Find avg jaw velocity for first condition trials with that delay
        jawvel.right{g} = medfilt1(tempjaw,10);         % Apply median filter

        gix = find(met.del_trialid{2}==g);              % Same thing for second condition
        tempjaw = nanmean(jaw_by_cond{2}(:,gix),2);
        jawvel.left{g} = medfilt1(tempjaw,10);
    end

    % Plot probability of jaw[Jaw] movement for each delay length
    colors = {[0 0 1],[1 0 0]};
    plotJawProb_HazardDel(taxis, jawvel,params)
    sesstitle = strcat('Prob of jaw movement for',{' '},anm,{' '},date);  % Name/title for session
    sgtitle(sesstitle,'FontSize',13)

    %%%% FIND CHOICE MODE %%%%
    rez.time = objs{1}.time;
    rez.condition = objs{1}.condition;
    rez.alignEvent = params.alignEvent;
    % Find CDearly (coding dimension during mid-sample)
    cond{1} = params.modecondition{1};
    cond{2} = params.modecondition{2};
    epoch = 'midsample';
    CD.early_mode = choiceMode(obj,met,cond,epoch,rez.alignEvent,'no');
    % Find CDlate (coding dimension during late delay period)
    epoch = 'latedelay';
    CD.late_mode = choiceMode(obj,met,cond,epoch,rez.alignEvent,'yes');
    clear cond
    % Find go mode
    cond{1} = params.modecondition{5};
    epoch = {'postgo','prego'};
    CD.go_mode = goMode(obj,met,cond,epoch,rez.alignEvent,'no');
    clear cond
    % Orthogonalize CD late to CD early
    CD = orthogModes(CD, obj);

    % Get the projection of specified conditions onto the choice mode
    smooth = 51;
    latent.early = getChoiceProj_byDelLength(CD.early_mode,smooth,delPSTH,params);
    latent.late = getChoiceProj_byDelLength(CD.late_mode,smooth,delPSTH,params);

    % Plot choice mode for each condition separated by delay length
    if strcmp(plotIndiv,'yes')
        figure();
        for e = 1:3                 % For CDearly, CDlate, and jaw prob...
            for d = 1:4             % For the first 4 delay lengths...
                if e==1             % CDearly
                    la = latent.early;
                    dur = 'CDearly';
                    s = d;
                    colors = {[0 0 1],[1 0 0]};
                elseif e==2         % CDlate
                    la = latent.late;
                    dur = 'CDlate';
                    colors = {[0 0 1],[1 0 0]};
                    s = d+4;
                elseif e==3         % Jaw stuff
                    la = jawvel;
                    s = d+8;
                    colors = {[0 0 0.9],[0.9 0 0]};
                end
                subplot(3,4,s)
                plot(taxis,la.left{d},'Color',colors{2},'LineWidth',2)
                hold on;
                plot(taxis,la.right{d},'Color',colors{1},'LineWidth',2)
                xlim([-1.4 2])
                xline(0,'LineStyle','--','LineWidth',1.3)
                xline(-1.3,'LineStyle',':','LineWidth',1.3)
                xline(params.delay(d),'LineStyle','-.','LineWidth',1.3)
                %legend('Left','Right','Delay onset','Sample onset','GoCue','Location','best')
                len = num2str(params.delay(d));
                if e==1 || e==2
                    subtitle = strcat(dur,';Delay length =',{' '},len);
                    ylabel('a.u.')
                elseif e==3
                    subtitle = strcat('Jaw prob; Delay length =',{' '},len);
                    ylabel('%')
                    xlabel('Time since delay onset')
                end
                title(subtitle)
            end
        end
        sesstitle = strcat({' '},anm,{' '},date);  % Name/title for session
        sgtitle(sesstitle,'FontSize',13)
    end
end
