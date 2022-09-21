clear,clc,close all

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
%% SET RUN PARAMS
params.alignEvent          = 'firstLick'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'
params.timeWarp = 0;
params.nLicks = 20;
params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val
params.moveThresh          = 0.15;      % What percentage of the delay period you want to use for identifying early move trials

% set conditions to calculate PSTHs for
% params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
% params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
% params.condition(1)     = {'L&hit&~stim.enable&~autowater&~early'};             % all 2AFC hits, no stim, aw off
% params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};              % all AW hits, no stim
params.condition(1)     = {'hit&~stim.enable&~autowater&~early'};             % all 2AFC hits, no stim, aw off
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};              % all AW hits, no stim


params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 51;

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
%meta = loadJEB15_ALMVideo(meta);

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

%% Combine cells from both probes in a session
anm = 'JEB15';
dates = {'2022-07-26','2022-07-27','2022-07-28'};       % Dates that you want to combine probe info from
%[meta, objs, params] = combineSessionProbes(meta,objs,params,anm,dates);
%% Adjust for older functions
for i = 1:numel(objs)
    meta(i).trialid = params.trialid{i};        % What used to be stored as "meta.trialid" is now "params.trialid"
    conditions = {1,2};

    for c = 1:numel(conditions)
        trix = meta(i).trialid{c};              % Get the condition-separated trial PSTHs from the trial data
        objs{i}.trialpsth_cond{c} = objs{i}.trialdat(:,:,trix);
    end
end
%%  Plot kinematic modes and kinematic measurements
params.kinfind = 'vel';
for gg = 4%1:length(meta)         % For all loaded sessions...
    clear orthproj orthModes proj kinmodes kinfns numTrials kin orthmode mode varexp
    
    obj = objs{gg};
    met = meta(gg);
    [anm,date,probenum,taxis] = getSessMeta(obj,met);
    
    % FIND KINEMATIC MODES %
    psthForProj = cat(3,obj.trialpsth_cond{1},obj.trialpsth_cond{2});       % Concatenated single trial PSTHs (all trials from cond 1 then all trials from cond 2)

    %%% GET FEATURE KINEMATICS (for all trial types) %%%
    conditions = {1,2};     
    kin = struct();
    kin.MEinterp = getME(obj,met,params,taxis,conditions);                  % Find interpolated motion energy
    kinfeat = findAllFeatKin(params, taxis, obj,conditions,met);            % Find all kinematic features (jaw, nose, tongue vel)
    kin.jawVel = kinfeat.jawVel; kin.noseVel = kinfeat.noseVel; kin.tongueVel = kinfeat.tongueVel;
    kin.tongueAngle = findTongueAngle(taxis, obj, met,params,conditions);   % Find tongue angle
   
    %%% DETERMINE WHICH TRIAL TYPE THE ANIMAL MOVES MORE ON %%%
    modecond = findMoreMove(conditions,kin,'jawVel',met,taxis);
    modekin = kinRestruct(kin,met,modecond);                                               % Only keep the kinematic features from the condition that you want (necessary for mode calc)
   
    %%% CALCULATE THE KINEMATIC MODES BASED ON THE SPECIFIED CONDITION %%%
    % Tongue angle is always calculated using both conditions %
    [mode,dat,proj,kinfns] = calcKinModes(modekin,obj,params,psthForProj,modecond,taxis);  % Time-points used to find modes are specified in the function

    % Orthogonalize modes
    orthmode = orthogonalizeModes(mode, kinfns);

    % Project data onto orthmodes
    for i = 1:numel(kinfns)
        orthproj.(kinfns{i}) = getProjection(psthForProj, orthmode.(kinfns{i}));
    end

    % VARIANCE EXPLAINED %
    varexp = findVarExplained(obj,orthmode,kinfns);
    
    %%% PLOT EVERYTHING %%%
    conditions = {1,2};
    plotModeFeat_SingleTrix(met,kinfns,proj,taxis,varexp,anm,date,probenum,kin,params)
    plotAvgModeProj(kinfns,conditions,met,proj,varexp,anm,date,probenum,taxis)

    disp('hi')
end
%%  Organize Variance Explained into correct structure for plotting
% [VEhalf,VEfull] = organizeVE(varexp_half, varexp_full);
% 
% %  Plot Variance explained across different ways of finding Kinematic modes
% plotVEScatterBar(VEhalf,VEfull)
%% FUNCTIONS

function MEinterp = getME(obj,met,params,taxis,conditions)
[met,mov,me] = assignEarlyTrials(obj,met,params);
[temp,~] = findInterpME(taxis,conditions, met,mov,me,params,obj);
MEinterp = [];
for c = 1:numel(conditions)
    MEinterp = [MEinterp, temp{c}];
end
end

function modekin = kinRestruct(kin,met,modecond)
kinfns = fieldnames(kin);
for f = 1:numel(kinfns)
    if modecond == 1
        ix = 1:length(met.trialid{1});
    else
        start = (length(met.trialid{1})+1);
        stop = start + met.trialid{2};
        ix = start:stop;
    end
    if ~strcmp(kinfns{f},'tongueAngle')
        modekin.(kinfns{f}) = kin.(kinfns{f})(:,ix);
    else
        modekin.(kinfns{f}) = kin.(kinfns{f});
    end
end
end

function varexp = findVarExplained(obj,orthmode,kinfns)
for i = 1:numel(kinfns)
        psth = obj.psth;
        datacov = cov([psth(:,:,1) ; psth(:,:,2)]);
        datacov(isnan(datacov)) = 0;
        eigsum = sum(eig(datacov));
        varexp.(kinfns{i}) = var_proj(orthmode.(kinfns{i}), datacov, eigsum);
end
end  % findVarExplained

function [VEhalf,VEfull] = organizeVE(varexp_half, varexp_full)
VEhalf.ME = [];  VEfull.ME = [];
VEhalf.jaw = [];  VEfull.jaw = [];
VEhalf.nose = [];  VEfull.nose = [];
for gg = 1:length(meta)
    VEhalf.ME = [VEhalf.ME, varexp_half(gg).MEinterp];  VEfull.ME = [VEfull.ME, varexp_full(gg).MEinterp];
    VEhalf.jaw = [VEhalf.jaw, varexp_half(gg).jawVel];  VEfull.jaw = [VEfull.jaw, varexp_full(gg).jawVel];
    VEhalf.nose = [VEhalf.nose, varexp_half(gg).noseVel];  VEfull.nose = [VEfull.nose, varexp_full(gg).noseVel];
end
end
%% PLOTTING FUNCTIONS
function plotModeFeat_SingleTrix(met,kinfns,orthproj,taxis,varexp,anm,date,probenum,kin,params)
numTrials = length(met.trialid{1})+length(met.trialid{2});
l1 = length(met.trialid{1});
l2 = l1+length(met.trialid{2});

for i=1:numel(kinfns)
    % Projections onto all modes
    %f(i) = subplot(numel(kinfns),2,(2*i-1));
    figure();
    subplot(1,2,1)
    imagesc(taxis,1:numTrials,orthproj.(kinfns{i})')
    colorbar()
    hold on;
    line([taxis(1),taxis(end)],[l1,l1],'Color','white','LineStyle','--')
    line([taxis(1),taxis(end)],[l2,l2],'Color','white','LineStyle','--')
    xlabel('Time since go-cue (s)')
    xlim([-2.45 2.4])
    VE = num2str(varexp.(kinfns{i}));
    plotitle = strcat(kinfns{i},'  Mode;',{' '},'VE = ',VE);
    title(plotitle)

    % Kinematic feature tracking
    %g(i) = subplot(numel(kinfns),2,(2*i));
    subplot(1,2,2)
    imagesc(taxis,1:numTrials,kin.(kinfns{i})')
    colorbar()
    if strcmp(params.kinfind,'vel')
        if strcmp(kinfns{i},'noseVel')
            caxis([0 1])
        elseif strcmp(kinfns{i},'tongueAngle')
            caxis([-2 2])
        else
            caxis([0 7])
        end
    else
        if strcmp(kinfns{i},'jawPos')
            caxis([0 20])  
        elseif strcmp(kinfns{i},'nosePos')
            caxis([0 10])  
        end
    end
    xlim([-2.45 2.4])
    hold on;
    line([taxis(1),taxis(end)],[l1,l1],'Color','white','LineStyle','--')
    line([taxis(1),taxis(end)],[l2,l2],'Color','white','LineStyle','--')
    xlabel('Time since go-cue (s)')
    title(kinfns{i},'  Feature Tracking');
end
sesstitle = strcat(anm,date,' ;  ','Probe ',probenum,'Kinematic Modes');  % Name/title for session
sgtitle(sesstitle,'FontSize',16)

end

function plotAvgModeProj(kinfns,conditions,met,orthproj,varexp,anm,date,probenum,taxis)

figure();
% colors = {[0 0 1],[1 0 0]};
% leg = {'Right','Left'};
colors = {'magenta','green'};
leg = {'L 2AFC','L AW'};
for k = 1:numel(kinfns)
    subplot(numel(kinfns),1,k)
    for c = 1:numel(conditions)
        cond = conditions{c};
        col = colors{c};
        if c==1
            maxmin = 1:length(met.trialid{c});
            nums = length(met.trialid{c});
        else
            start = nums(c-1)+1;
            fin = start+length(met.trialid{c})-1;
            maxmin = start:fin;
            nums = [nums,fin];
        end
        projtouse = orthproj.(kinfns{k})(:,maxmin);
        avgproj = mean(projtouse,2);
        plot(taxis,avgproj,'Color',col,'LineWidth',2)
        hold on;
    end
    VE = num2str(varexp.(kinfns{k}));
    plotitle = strcat(kinfns{k},'  Mode;',{' '},'VE = ',VE);
    legend(leg)
    title(plotitle)
end
sesstitle = strcat(anm,date,' ;  ','Probe ',probenum,'Kinematic Modes');  % Name/title for session
sgtitle(sesstitle,'FontSize',16)
end

function plotVEScatterBar(VEhalf,VEfull)
figure();
x = [1,2,4,5,7,8];
y = [mean(VEhalf.ME),mean(VEfull.ME,'omitnan'),mean(VEhalf.jaw),mean(VEfull.jaw,'omitnan'),mean(VEhalf.nose),mean(VEfull.nose,'omitnan')];
b= bar(x,y);
b.FaceColor = [0.75 0.75 0.75]; hold on;

scatter(1,VEhalf.ME,'yellow','filled'); scatter(2,VEfull.ME,'black','filled')
scatter(4,VEhalf.jaw,'yellow','filled');  scatter(5,VEfull.jaw,'black','filled')
scatter(7,VEhalf.nose,'yellow','filled');  scatter(8,VEfull.nose,'black','filled')
yy = [VEhalf.ME',VEfull.ME']; plot(x(1:2),yy(:,1:2),'Color','black')
yy = [VEhalf.jaw',VEfull.jaw']; plot(x(3:4),yy(:,1:2),'Color','black')
yy = [VEhalf.nose',VEfull.nose']; plot(x(5:6),yy(:,1:2),'Color','black')

xlim([0 9])
ylabel('Variance explained (out of 1)')
title('Difference in VE using 1/2 trials vs all trials')
end