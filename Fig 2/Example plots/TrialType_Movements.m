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
params.condition(1)     = {'hit&~stim.enable&~autowater&~early'};             % all 2AFC hits, no stim, aw off
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};              % all AW hits, no stim


params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 51;

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
for gg = 2 %:length(meta)         % For all loaded sessions...
    clear orthproj orthModes proj kinmodes kinfns numTrials kin orthmode mode varexp
    
    obj = objs{gg};
    met = meta(gg);
    [anm,date,probenum,taxis] = getSessMeta(obj,met);
    
    psthForProj = cat(3,obj.trialpsth_cond{1},obj.trialpsth_cond{2});       % Concatenated single trial PSTHs (all trials from cond 1 then all trials from cond 2)

    %%% GET FEATURE KINEMATICS (for all trial types) %%%
    conditions = {1,2};     
    kin = struct();
    kinfeat = findAllFeatKin(params, taxis, obj,conditions,met);            % Find all kinematic features (jaw, nose, tongue vel)

    %plotModeFeat_SingleTrix(met,kinfeat,taxis,anm,date,probenum,params)
    plotAvgModeProj(conditions,met,kinfeat,anm,date,probenum,taxis)
end
%% PLOTTING FUNCTIONS
function plotModeFeat_SingleTrix(met,kinfeat,taxis,anm,date,probenum,params)
numTrials = length(met.trialid{1})+length(met.trialid{2});
l1 = length(met.trialid{1});
l2 = l1+length(met.trialid{2});
kinfns = fieldnames(kinfeat);
for i=1:length(kinfns)
    % Kinematic feature tracking
    figure();
    imagesc(taxis,1:numTrials,kinfeat.(kinfns{i})')
    colorbar()
    if strcmp(params.kinfind,'vel')
        if strcmp(kinfns{i},'noseVel')
            caxis([0 1])
        elseif strcmp(kinfns{i},'tongueAngle')
            caxis([-2 2])
        else
            caxis([0 6.5])
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

function plotAvgModeProj(conditions,met,kinfeat,anm,date,probenum,taxis)
figure();
kinfns = fieldnames(kinfeat);
colors = {[0 0 1],[1 0 0]};
leg = {'Right','Left'};
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
        kintouse = kinfeat.(kinfns{k})(:,maxmin);
        avgproj = mean(kintouse,2,'omitnan');
        plot(taxis,avgproj,'Color',col,'LineWidth',2)
        hold on;
    end
    plotitle = strcat(kinfns{k});
    legend(leg)
    title(plotitle)
    xlim([-2 2])
end
sesstitle = strcat(anm,date,' ;  ','Probe ',probenum,'Kinematic Modes');  % Name/title for session
sgtitle(sesstitle,'FontSize',16)
end