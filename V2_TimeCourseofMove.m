clear,clc,close all

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'
params.timeWarp = 0;
params.nLicks = 20;
params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val
params.context             = '2AFC';
params.moveThresh          = 0.15;      % What percentage of the delay period you want to use for identifying early move trials

% set conditions to calculate PSTHs for
params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
% params.condition(1)     = {'hit&~stim.enable&~autowater&~early'};             % all 2AFC hits, no stim, aw off
% params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};              % all AW hits, no stim


params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 31;

% Params for finding kinematic modes
params.fcut = 50;          % smoothing cutoff frequency
params.cond = 1:2;         % which conditions to use to find mode
params.method = 'xcorr';   % 'xcorr' or 'regress' (basically the same)
params.fa = false;         % if true, reduces neural dimensions to 10 with factor analysis

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
meta = loadJEB14_ALMVideo(meta);
meta = loadJEB15_ALMVideo(meta);

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
%%  Plot kinematic modes and kinematic measurements
params.kinfind = 'vel';
rxntime = cell(1,length(meta));
for gg = 1:length(meta)         % For all loaded sessions...
    clear orthproj orthModes proj kinmodes kinfns numTrials orthmode mode varexp

    obj = objs{gg};
    met = meta(gg);
    [anm,date,probenum,taxis] = getSessMeta(obj,met);

    %%% GET FEATURE KINEMATICS %%%
    conditions = {1,2};
    kin(gg).MEinterp = getME(obj,met,params,taxis,conditions);                  % Find interpolated motion energy
    kinfeat = findAllFeatKin(params, taxis, obj,conditions,met);
    kin(gg).jawVel = kinfeat.jawVel; kin(gg).noseVel = kinfeat.noseVel; kin(gg).tongueVel = kinfeat.tongueVel;
    kin(gg).tongueAngle = findTongueAngle(taxis, obj, met,params,conditions);
    
    %%% GET REACTION TIME ON EACH TRIAL %%%
    rxntime{gg} = getRxnTimes(met,conditions,obj);
end
%%  Find average jaw movement on each trial during the whole delay period
e1 = find(taxis>-0.85,1,'first');
e2 = find(taxis>-0.05,1,'first');
moveDelay = findAvgMove(meta,kin,e1,e2);
%%
%%%% Find avg movement during first third of session, second third, and
%%%% last third
numDivisions = 4;
for gg = 1:length(meta)
    numR = length(moveDelay.right{gg});  numL = length(moveDelay.left{gg});
    divideR = round(numR/numDivisions); divideL = round(numL/numDivisions);
    cntR = 1; cntL = 1;
    for d = 1:numDivisions
        if cntR==1
            averageR{d}(gg) = mean(moveDelay.right{gg}(1:divideR),2,'omitnan');  averageL{d}(gg) = mean(moveDelay.left{gg}(1:divideL),2,'omitnan');
            cntR = divideR+1; cntL = divideL+1;
        else
            if cntR+divideR>numR
                stopixR = numR;
            else
                stopixR = cntR+divideR;
            end
            if cntL+divideL>numL
                stopixL = numL;
            else
                stopixL = cntL+divideL;
            end
            averageR{d}(gg) = mean(moveDelay.right{gg}(cntR:stopixR),2,'omitnan');  averageL{d}(gg) = mean(moveDelay.left{gg}(cntL:stopixL),2,'omitnan');
            cntR = cntR+divideR+1; cntL = cntL+divideL+1;
           
        end
    end
end
%%
plotAvgKin_segments(averageL,averageR,numDivisions)
%%
clear rt
figure();
md = [];
rt = [];
for gg = 1:length(meta)
    move = [moveDelay.right{gg},moveDelay.left{gg}];
    md = [md,move];
    rt = [rt,rxntime{gg}];
end
scatter(md,rt)
xlabel("Avg ME during delay")
ylabel('Reaction time (s)')
%% FUNCTIONS
function moveDelay = findAvgMove(meta,kin,e1,e2)
moveDelay.right = cell(1,length(meta));
moveDelay.left = cell(1,length(meta));
for gg = 1:length(meta)
    met = meta(gg);
    rightix = 1:length(met.trialid{1});
    leftstart = (length(met.trialid{1})+1);
    tempR = kin(gg).MEinterp(e1:e2,rightix);  tempR = mean(tempR,1); moveDelay.right{gg} = tempR;
    tempL = kin(gg).MEinterp(e1:e2,leftstart:end); tempL = mean(tempL,1); moveDelay.left{gg} = tempL;
end
end

function MEinterp = getME(obj,met,params,taxis,conditions)
[met,mov,me] = assignEarlyTrials(obj,met,params);
[temp,~] = findInterpME(taxis,conditions, met,mov,me,params,obj);
MEinterp = [];
for c = 1:numel(conditions)
    MEinterp = [MEinterp, temp{c}];
end
end

function rxntime = getRxnTimes(met,conditions,obj)
nTrials = length(met.trialid{1}) + length(met.trialid{2});
rxntime = NaN(1,nTrials);
cnt = 0;
for c = 1:numel(conditions)
    cond = conditions{c};
    trix = met.trialid{cond};
    for t = 1:length(trix)
        cnt = cnt+1;
        currtrial = trix(t);
        allLicks = [obj.bp.ev.lickL{currtrial}, obj.bp.ev.lickR{currtrial}];
        temp = allLicks>obj.bp.ev.goCue(currtrial); allLicks = allLicks(temp);
        allLicks = sort(allLicks,'ascend');
        firstLick = allLicks(1);
        rt = firstLick - obj.bp.ev.goCue(currtrial);
        rxntime(cnt) = rt;
    end
end
end

function plotAvgKin_segments(averageL,averageR,numDivisions)
figure();
x = [1,2,3,4,6,7,8,9];
yL = []; yR = [];
for d = 1:numDivisions
    yL = [yL, mean(averageL{d})]; yR = [yR, mean(averageR{d})];
end
y = [yL, yR];
b= bar(x,y);
b.FaceColor = [0.75 0.75 0.75]; hold on;

scatter(1,averageL{1},'red','filled'); scatter(2,averageL{2},'red','filled')  ; scatter(3,averageL{3},'red','filled'); scatter(4,averageL{4},'red','filled')
scatter(6,averageR{1},'blue','filled'); scatter(7,averageR{2},'blue','filled')  ; scatter(8,averageR{3},'blue','filled'); scatter(9,averageR{4},'blue','filled')
yy = [averageL{1}',averageL{2}',averageL{3}',averageL{4}']; plot(x(1:4),yy(:,1:4),'Color','black')
yy = [averageR{1}',averageR{2}',averageR{3}',averageR{4}']; plot(x(5:8),yy(:,1:4),'Color','black')


xlim([0 10])
ylabel('Average motion energy')
title('ME throughout session')
end