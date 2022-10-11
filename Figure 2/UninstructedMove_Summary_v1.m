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
anm = 'JEB15';
dates = {'2022-07-26','2022-07-27','2022-07-28'};       % Dates that you want to combine probe info from
[meta, objs, params] = combineSessionProbes(meta,objs,params,anm,dates);
%% EXAMPLE HEATMAP OF JAW VELOCITY ON SINGLE TRIALS--SEPARATED BY TRIAL TYPE

sesh = 13;           % Maybe 1,4
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

% Sort the trials by jaw velocity during the late delay period
% Find the average jaw velocity during specified time points (on each trial)
    startix = find(taxis>=-0.4, 1, 'first');
    stopix = find(taxis<=-0.05, 1, 'last');
    val.right = nanmean(jaw_by_cond{1}(startix:stopix, :), 1);
    val.left = nanmean(jaw_by_cond{2}(startix:stopix, :), 1);
    nanix = find(isnan(val.right)); nanix = find(isnan(val.left));                   % Get rid of trials where jaw velocity is always NaN
    val.right(nanix) = [];  val.left(nanix) = [];
    
    % Sort the average jaw velocities in descending order and save the trial
    % order
    sort_by_cond = cell(1,numel(conditions));
    [sort_by_cond{1}, six1] = sort(val.right, 'descend'); jaw_by_cond{1} = jaw_by_cond{1}(:,six1);
    [sort_by_cond{2}, six2] = sort(val.left, 'descend'); jaw_by_cond{2} = jaw_by_cond{2}(:,six2);

% Plot
figure();
numTrixPlot = 30;
rangetoPlot = 1:numTrixPlot;
tGo = find(taxis==0);
for i=1:numel(conditions)
    if i==1
        subplot(1,2,2);
    elseif i==2
        subplot(1,2,1);
    end
    imagesc(taxis,1:numTrixPlot,mySmooth(jaw_by_cond{i}(:,1:numTrixPlot)',5)); caxis([0 2]); colormap("hot");
    go = 0; 
    delstart = -0.9;
    sampstart = delstart-1.3;
    line([go,go],[0,numTrixPlot+0.5],'Color','white','LineStyle','--')
    line([delstart,delstart],[0,numTrixPlot+0.5],'Color','white','LineStyle','--')
    line([sampstart,sampstart],[0,numTrixPlot+0.5],'Color','white','LineStyle','--')
    xlim([taxis(10), 0])
    xlabel('Time before go-cue (s)','FontSize',13)
    ylabel('Trials','FontSize',13)
    if i==1
        title('Right trials')
    elseif i==2
        title('Left trials')
    end
    c=colorbar;
    if strcmp(params.jawMeasure,'sideJaw')
        ylabel(c,'Jaw velocity','FontSize',12,'Rotation',90);
    else
        ylabel(c,'Motion Energy','FontSize',12,'Rotation',90);
    end
end
figtitle =  strcat('Example Trials from',anm,date,' ;  ','Probe ',probenum);  % Name/title for session
sgtitle(figtitle,'FontSize',16)
%% EXAMPLE SESSION OF PROBABILITY OF JAW MOVEMENT--SEPARATED BY TRIAL TYPE
obj = objs{sesh};     
met = meta(sesh);

anm = obj.pth.anm;                  % Animal name
date = obj.pth.dt;                  % Session date
probenum = string(met.probe);       % Which probe was used

conditions = {1,2};
colors = {[0 0 1],[1 0 0]};

% Plot
figure();
if strcmp(params.jawMeasure,'sideJaw')
    plotJawProb_SessAvg(obj,met,conditions,colors,taxis,'no',params)
elseif strcmp(params.jawMeasure,'MotionEnergy')
    plot(taxis,mySmooth(mean(jaw_by_cond{1},2),31),'Color',colors{1},'LineWidth',2.5); hold on;
    plot(taxis,mySmooth(mean(jaw_by_cond{2},2),31),'Color',colors{2},'LineWidth',2.5);
    xlabel('Time (s) before go-cue')
    ylabel('Motion Energy')
end
legend('Right','Left','Location','best')

% Add lines at trial times
go = 0;
trix = met.trialid{1}(1);
del = obj.bp.ev.goCue(trix)-obj.bp.ev.delay(trix);
delstart = 0-del;
sampstart = delstart-1.3;
xline(go,'Color','black','LineStyle','--')
xline(delstart,'Color','black','LineStyle','--')
xline(sampstart,'Color','black','LineStyle','--')
xlim([-2.3 0])

figtitle =  strcat('Example w/ selectivity',anm,date,' ;  ','Probe ',probenum);  % Name/title for session
title(figtitle,'FontSize',16)
%% EXAMPLE SESSION OF CDLate projections, SEPARATED BY TRIAL TYPE
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

% Plot
figure();
colors = {[0 0.4470 0.7410],[0.6350 0.0780 0.1840]};
for c = 1:numel(cond)
    temp = cond(c);
    col = colors{c};
    subplot(2,1,1)
    plot(taxis,mySmooth(rez(sesh).cdLate_latent(:,temp),51),'Color',col,'LineWidth',3)
    hold on;
    set(gca, 'YDir','reverse')
    ylabel('a.u.')
    col = 'black';
    addTrialLines(col,met,obj)
    xlim([-2.3 0])
    xlabel('Time since go-cue')
    
    subplot(2,1,2)
    plot(taxis,mySmooth(mean(jaw_by_cond{temp},2),51),'Color',colors{temp},'LineWidth',2.5); hold on;
    ylabel('M.E.')
    addTrialLines(col,met,obj)
    xlim([-2.3 0])
    xlabel('Time since go-cue')
end


legend('Right','Left','Location','best')

figtitle =  strcat('Example CDchoice',anm,date,' ;  ','Probe ',probenum);  % Name/title for session
title(figtitle,'FontSize',16)
%% EXAMPLE SESSION--SCATTER PLOT OF SINGLE TRIAL JAW VELOCITY VS CHOICE CD

obj = objs{sesh};     % 11th data object = JEB7, 04-29 (Classic sesh)
met = meta(sesh);

anm = obj.pth.anm;                  % Animal name
date = obj.pth.dt;                  % Session date
probenum = string(met.probe);       % Which probe was used

rez(sesh).time = objs{1}.time;
rez(sesh).alignEvent = params.alignEvent;

% Project single trials onto choice mode
cd = rez(sesh).cdLate_mode;
latent = getTrialLatents(obj,cd,conditions,met);
lat_choice = [];
jaw = [];
for c = 1:numel(conditions)
    lat_choice = [lat_choice,latent{c}];
    jaw = [jaw,jaw_by_cond{c}];
end

% Define time intervals: Time frame for late delay period(from -0.4 before go-cue to -0.1)
late_start = find(taxis>=-0.4, 1, 'first');
late_stop = find(taxis<=-0.05, 1, 'last');
lateDelay = late_start:late_stop;

% Get jaw velocity and activity mode averages for late delay
timeInt = lateDelay;
jawVel_late = getAverages(timeInt,jaw);
Choice_late = getAverages(timeInt,lat_choice);

% Make scatter plot
conditions = 1:2;               % Look only at correct left and right hits during 2AFC
figure();
colors = {[0 0 0],[0 0 0]};
ActivityMode_Jaw_Scatter(jawVel_late,Choice_late,conditions,met,colors,obj,params);

nanix = find(isnan(jawVel_late));                    % Find indices where the jaw vel is a NaN
jv = jawVel_late(~isnan(jawVel_late));               % Get rid of the NaN values
Choice_late(nanix) = [];                             % Indices that were a NaN for jaw vel, get rid of those indices in the choice mode as well
ch = Choice_late;
R = corr2(jv,ch);                           % Calculate the correlation coefficient between these two variables
R = num2str(R);
coeff = polyfit(jv,ch,1);                   % Find the line of best fit
hline = refline(coeff);
hline.LineStyle = '--';
hline.Color = 'k';
str = strcat('R^2 =',R);
lgd = legend('Right','Left',str);
lgd.FontSize = 11;
lgd.Location = 'best';
if strcmp(params.jawMeasure,'sideJaw')
    xlabel('Avg Jaw Velocity','fontsize',14)
else
    xlabel('Avg Motion Energy','fontsize',14)
end
%xlim([0 1.2])
ylabel('Avg choice mode','fontsize',14)
figtitle =  strcat('Example Single trial correlations',anm,date,' ;  ','Probe ',probenum);  % Name/title for session
sgtitle(figtitle)

disp('hi');