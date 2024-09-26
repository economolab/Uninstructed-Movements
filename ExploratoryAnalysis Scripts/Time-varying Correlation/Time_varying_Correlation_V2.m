clear,clc,close all

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
%% SET RUN PARAMS
params.alignEvent          = 'firstLick'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'
params.timeWarp = 0;
params.nLicks = 20;
params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val
params.context             = '2AFC';

% set conditions to calculate PSTHs for
params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
% params.condition(1)     = {'hit&~stim.enable&~autowater&~early'};         % all 2AFC hits, no stim, aw off
% params.condition(1) = {'hit&~stim.enable&autowater&~early'};              % all AW hits, no stim


params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;
params.moveThresh = 0.15;

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
meta = loadJEB15_ALMVideo(meta);
meta = loadJEB14_ALMVideo(meta);
% meta = loadJEB14_M1TJVideo(meta);

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
inclJEB15 = 'yes';
if strcmp(inclJEB15,'yes')
    anm = 'JEB15';
    dates = {'2022-07-26','2022-07-27','2022-07-28'};       % Dates that you want to combine probe info from
    [meta, objs, params] = combineSessionProbes(meta,objs,params,anm,dates);
end
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

for mm = 1:length(meta)         % For all loaded sessions...
    clear orthproj orthModes proj kinmodes kinfns numTrials kin orthmode mode varexp
    %gg = touse(mm);
    gg =mm;

    obj = objs{gg};
    met = meta(gg);
    [anm,date,probenum,taxis] = getSessMeta(obj,met);

    psthForProj = cat(3,obj.trialpsth_cond{1},obj.trialpsth_cond{2});       % Concatenated single trial PSTHs (all trials from cond 1 then all trials from cond 2)
    %psthForProj = cat(3,obj.trialpsth_cond{1});

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

    % Find time-varying correlation for a kinematic feature %
    kinfeat = 'jawVel';
    window = 0.05;                                                              % Time window during which to find the correlation between FR and kin feature                                                         % Sort according to mean correlation prior to response per
    e1 = find(taxis>-0.9,1,'first');
    e2 = find(taxis>-0.1,1,'first');
        sortrange = 1:e2;
        XTrialtscorr = findXTrialTSCorr(psthForProj, kin.(kinfeat),taxis,window,sortrange);
    
        plotTSCorr(XTrialtscorr,kinfeat,taxis,params)
        sesstitle = strcat(anm,date);
        title(sesstitle)

    % Find Jaw Modes across each time point for the current session
    [blah(gg).tsModes,plusmin] = findTSModes(psthForProj, kin.(kinfeat),taxis,window);
end
%% Find a correlation matrix for all of the Jaw Modes found across time
for gg = 1:length(meta)
    modes = blah(gg).tsModes;           
    corrmatrix = findCorrMatrix(modes);
    blah(gg).corrmatrix = corrmatrix;
end
%%
% Plot correlation matrix for each session
for gg = 1:length(meta)
    corrmatrix = blah(gg).corrmatrix;
    obj = objs{gg};
    met = meta(gg);
    [anm,date,probenum,taxis] = getSessMeta(obj,met);
    plotCorrMatrix(taxis,plusmin,corrmatrix,params)
    sesstitle = strcat(anm,date);
    title(sesstitle)
end
%% Compare the contributions of all neurons to Jaw modes found at different time-points
% Specify time points for a given epoch within the trial
times.lateDel = e1:e2;      % Late delay
times.Resp = 700:790;       % Response period/consumption
times.MoveInit = 500:590;   % Movement initiation

% For every session...
for gg = 1:length(meta)
    lateDelWeights = mean((blah(gg).tsModes(:,times.lateDel)),2);        % Get the average weight for each neuron in the jaw modes found during the late delay period
    RespWeights = mean((blah(gg).tsModes(:,times.Resp)),2);              % Late response period
    MoveInitWeights = mean((blah(gg).tsModes(:,times.MoveInit)),2);      % During movement initiation/very early response period
    
    % Find how orthogonal each of these average weights are
    % TO-DO: ask the best way to do this 
    dotpro(gg).late_move = dot(lateDelWeights,MoveInitWeights);
    dotpro(gg).late_resp = dot(lateDelWeights,RespWeights);
    dotpro(gg).move_resp = dot(MoveInitWeights,RespWeights);
end
%% Find projections onto jaw modes found at different points throughout the trials
for gg = 1:length(meta)
    obj = objs{gg};

    latents(gg).lateDel = NaN(length(taxis),numel(conditions));
    latents(gg).MoveInit = NaN(length(taxis),numel(conditions));
    latents(gg).Resp = NaN(length(taxis),numel(conditions));

    temp.lateDel = NaN(length(taxis),length(times.lateDel),numel(conditions)); % Time in trial x Num time points to be averaged over x num conditions
    temp.MoveInit = NaN(length(taxis),length(times.MoveInit),numel(conditions));
    temp.Resp = NaN(length(taxis),length(times.Resp),numel(conditions));
    
    % For each point in time within a given epoch, find the projection onto
    % the jaw mode for that time point
    for c = 1:numel(conditions)
        for t = 1:length(times.lateDel)
            temp.lateDel(:,t,c) = obj.psth(:,:,c)*blah(gg).tsModes(:,times.lateDel(t));
        end

        for t = 1:length(times.MoveInit)
            temp.MoveInit(:,t,c) = obj.psth(:,:,c)*blah(gg).tsModes(:,times.MoveInit(t));
            temp.Resp(:,t,c) = obj.psth(:,:,c)*blah(gg).tsModes(:,times.Resp(t));
        end
        % Find the average latent for that time epoch (average across all
        % of the latents that you found)
        latents(gg).lateDel(:,c) = mean(temp.lateDel(:,:,c),2,'omitnan');
        latents(gg).MoveInit(:,c) = mean(temp.MoveInit(:,:,c),2,'omitnan');
        latents(gg).Resp(:,c) = mean(temp.Resp(:,:,c),2,'omitnan');
    end
end
%% Plot the average latents for each time epoch (for each session)
colors = {[0 0.4470 0.7410],[0.6350 0.0780 0.1840]};
for gg = 1:length(meta)
    sesstitle = strcat(meta(gg).anm,meta(gg).date);
    figure();
    for c = 1:numel(conditions)
        subplot(3,1,1)
        plot(taxis,latents(gg).lateDel(:,c),'Color',colors{c},'LineWidth',2); hold on
        title('Late Del')

        subplot(3,1,2)
        plot(taxis,latents(gg).MoveInit(:,c),'Color',colors{c},'LineWidth',2); hold on
        title('Move Init')

        subplot(3,1,3)
        plot(taxis,latents(gg).Resp(:,c),'Color',colors{c},'LineWidth',2); hold on
        title('Response')
    end
    sgtitle(sesstitle)
    
end
%% Find the average latents + confidence intervals across all sessions
[cdlateDel, cdMoveInit, cdResp] = AvgKinModes_AcrossSess(latents,taxis,meta);
[upperci,lowerci] = findConfInt_KM(meta,taxis,conditions,cdlateDel,cdMoveInit,cdResp);
%% Plot the average latents across all sessions
colors = {[0 0.4470 0.7410],[0.6350 0.0780 0.1840]};
figure();
for c = 1:numel(conditions)
    subplot(3,1,1)
    plot(taxis,cdlateDel.avg(:,c),'Color',colors{c},'LineWidth',2); hold on
    patch([taxis(10:end) fliplr(taxis(10:end))],[lowerci.lateDel(9:end,c)' fliplr(upperci.lateDel(9:end,c)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
    title('Late Del')

    subplot(3,1,2)
    plot(taxis,cdMoveInit.avg(:,c),'Color',colors{c},'LineWidth',2); hold on
    patch([taxis(10:end) fliplr(taxis(10:end))],[lowerci.MoveInit(9:end,c)' fliplr(upperci.MoveInit(9:end,c)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
    title('Move Init')

    subplot(3,1,3)
    plot(taxis,cdResp.avg(:,c),'Color',colors{c},'LineWidth',2); hold on
    patch([taxis(10:end) fliplr(taxis(10:end))],[lowerci.Resp(9:end,c)' fliplr(upperci.Resp(9:end,c)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
    title('Response')
end
%% FUNCTIONS
function [upperci,lowerci] = findConfInt_KM(meta,taxis,conditions,cdlateDel,cdMoveInit,cdResp)
nSessions = length(meta);
upperci.lateDel = NaN(length(taxis)-1,2); lowerci.lateDel = NaN(length(taxis)-1,2);
upperci.MoveInit = NaN(length(taxis)-1,2); lowerci.MoveInit = NaN(length(taxis)-1,2);
upperci.Resp = NaN(length(taxis)-1,2); lowerci.Resp = NaN(length(taxis)-1,2);
for c = 1:numel(conditions)
    upperci.lateDel(:,c) = cdlateDel.avg(2:end,c)+1.96*(cdlateDel.std(2:end,c)/nSessions);  % Find the upper 95% confidence interval for each condition
    lowerci.lateDel(:,c) = cdlateDel.avg(2:end,c)-1.96*(cdlateDel.std(2:end,c)/nSessions);  % Find lower 95% condifence interval for each condition

    upperci.MoveInit(:,c) = cdMoveInit.avg(2:end,c)+1.96*(cdMoveInit.std(2:end,c)/nSessions);  
    lowerci.MoveInit(:,c) = cdMoveInit.avg(2:end,c)-1.96*(cdMoveInit.std(2:end,c)/nSessions);  

    upperci.Resp(:,c) = cdResp.avg(2:end,c)+1.96*(cdResp.std(2:end,c)/nSessions);  
    lowerci.Resp(:,c) = cdResp.avg(2:end,c)-1.96*(cdResp.std(2:end,c)/nSessions);   
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
%% PLOTTING FUNCTIONS
function plotTSCorr(tscorr,kinfeat,taxis,params)
figure();
imagesc(taxis(7:end-3),1:size(tscorr,1),tscorr);
colorbar
colormap('jet')
xlab = strcat('Time since','  ',params.alignEvent,'  (s)');
xlabel(xlab)
ylabel('Cells')
plottitle = strcat('Time-varying correlation with',kinfeat);
title(plottitle)
end  % plotTSCorr

function plotCorrMatrix(taxis,plusmin,corrmatrix,params)
figure();
imagesc(taxis(plusmin:end-plusmin),taxis(plusmin:end-plusmin),corrmatrix);
c = colorbar; ylabel(c,'Correlation between Jaw Modes','FontSize',13,'Rotation',270); c.Label.Position(1) = 4;
colormap('jet')
xline(0,'LineStyle','--','Color','black','LineWidth',2)
xline(-0.9,'LineStyle','--','Color','black','LineWidth',2)
xline(-2.4,'LineStyle','--','Color','black','LineWidth',2)
yline(0,'LineStyle','--','Color','black','LineWidth',2)
yline(-0.9,'LineStyle','--','Color','black','LineWidth',2)
yline(-2.4,'LineStyle','--','Color','black','LineWidth',2)
xlab = strcat('Time since','  ',params.alignEvent,'(s)');
xlabel(xlab,'FontSize',15)
ylabel(xlab,'FontSize',15)
end

