clear,clc,close all

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'
params.timeWarp = 0;
params.nLicks = 20;
params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val
params.context             = '2AFC';

% set conditions to calculate PSTHs for
% params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
% params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
params.condition(1)     = {'hit&~stim.enable&~autowater&~early'};             % all 2AFC hits, no stim, aw off
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};              % all AW hits, no stim


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
%%  Plot kinematic modes and kinematic measurements
params.kinfind = 'vel';
for gg = 1:length(meta)         % For all loaded sessions...
    clear orthproj orthModes proj kinmodes kinfns numTrials kin orthmode mode varexp
    
    obj = objs{gg};
    met = meta(gg);
    [anm,date,probenum,taxis] = getSessMeta(obj,met);

    %%% GET FEATURE KINEMATICS %%%
    kin = struct();
    kin = findAllFeatKin(params, taxis, obj,conditions,met);
    kin.tongueAngle = findTongueAngle(taxis, obj, met,params,conditions);

    % FIND KINEMATIC MODES %
    psthForProj = cat(3,obj.trialpsth_cond{1},obj.trialpsth_cond{2});       % Concatenated single trial PSTHs (all trials from cond 1 then all trials from cond 2)
    
    e1 = find(taxis>-0.5,1,'first');
    e2 = find(taxis>-0.05,1,'first');
    params.tix = e1:e2;        % time points to use when finding mode (LATE DELAY)
    [mode,dat,proj,kinfns] = calcKinModes(kin,obj,params,psthForProj);


    % Orthogonalize modes
    orthmode = orthogonalizeModes(mode, kinfns);

    % Project data onto orthmodes
    for i = 1:numel(kinfns)
        orthproj.(kinfns{i}) = getProjection(psthForProj, orthmode.(kinfns{i}));
    end

    % VARIANCE EXPLAINED %
    varexp = findVarExplained(obj,orthmode,kinfns);
    
    %%% PLOT EVERYTHING %%%
    plotModeFeat_SingleTrix(met,kinfns,orthproj,taxis,varexp,anm,date,probenum,kin,params)
    plotAvgModeProj(kinfns,conditions,met,orthproj,varexp,anm,date,probenum,taxis,params)

end
%% FUNCTIONS
function varexp = findVarExplained(obj,orthmode,kinfns)
for i = 1:numel(kinfns)
        psth = obj.psth;
        datacov = cov([psth(:,:,1) ; psth(:,:,2)]);
        datacov(isnan(datacov)) = 0;
        eigsum = sum(eig(datacov));
        varexp.(kinfns{i}) = var_proj(orthmode.(kinfns{i}), datacov, eigsum);
end
end  % findVarExplained
%% PLOTTING FUNCTIONS
function plotModeFeat_SingleTrix(met,kinfns,orthproj,taxis,varexp,anm,date,probenum,kin,params)
numTrials = length(met.trialid{1})+length(met.trialid{2});
l1 = length(met.trialid{1});
l2 = l1+length(met.trialid{2});

ff = figure();
ff.WindowState = 'maximized';
for i=1:numel(kinfns)
    % Projections onto all modes
    f(i) = subplot(numel(kinfns),2,(2*i-1));
    imagesc(taxis,1:numTrials,orthproj.(kinfns{i})')
    colorbar(f(i))
    hold on;
    line([taxis(1),taxis(end)],[l1,l1],'Color','white','LineStyle','--')
    line([taxis(1),taxis(end)],[l2,l2],'Color','white','LineStyle','--')
    xlabel('Time since go-cue (s)')
    xlim([-2.45 2.4])
    VE = num2str(varexp.(kinfns{i}));
    plotitle = strcat(kinfns{i},'  Mode;',{' '},'VE = ',VE);
    title(plotitle)

    % Kinematic feature tracking
    g(i) = subplot(numel(kinfns),2,(2*i));
    imagesc(taxis,1:numTrials,kin.(kinfns{i})')
    colorbar(g(i))
    if strcmp(params.kinfind,'vel')
        if strcmp(kinfns{i},'noseVel')
            caxis([0 1])
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

function plotAvgModeProj(kinfns,conditions,met,orthproj,varexp,anm,date,probenum,taxis,params)
figure();
if strcmp(params.context,'AW')
    colors = {[0.89 0.47 0.2],'green'}; 
    leg = {'2AFC','AW'};
else
    colors = {[0 0 1],[1 0 0]};
    leg = {'Right','Left'};
end
for k = 1:numel(kinfns)
    subplot(numel(kinfns),1,k)
    for c = 1:numel(conditions)
        cond = conditions{c};
        col = colors{c};
        if c==1
            range = 1:length(met.trialid{c});
            nums = length(met.trialid{c});
        else
            start = nums(c-1)+1;
            fin = start+length(met.trialid{c})-1;
            range = start:fin;
            nums = [nums,fin];
        end
        projtouse = orthproj.(kinfns{k})(:,range);
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