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

    psthForProj = cat(3,obj.trialpsth_cond{1},obj.trialpsth_cond{2});           % Concatenated single trial PSTHs (all trials from cond 1 then all trials from cond 2)

    % Find time-varying correlation for a kinematic feature %
    kinfeat = 'jawVel';
    window = 0.05;                                                              % Time window during which to find the correlation between FR and kin feature                                                         % Sort according to mean correlation prior to response per
    e1 = find(taxis>-0.5,1,'first');
    e2 = find(taxis>-0.05,1,'first');
    %     sortrange = 1:e2;
    %     XTrialtscorr = findXTrialTSCorr(psthForProj, kin.(kinfeat),taxis,window,sortrange);
    %
    %     plotTSCorr(XTrialtscorr,kinfeat,taxis)

    % Find Jaw Modes across each time point for the current session
    [blah(gg).tsModes,plusmin] = findTSModes(psthForProj, kin.(kinfeat),taxis,window);
end
%% Find a correlation matrix for all of the Jaw Modes found across time
for gg = 1:length(meta)
    modes = blah(gg).tsModes;           
    corrmatrix = findCorrMatrix(modes);
    blah(gg).corrmatrix = corrmatrix;
end
%% Plot correlation matrix
for gg = 1:length(meta)
    corrmatrix = blah(gg).corrmatrix;
    plotCorrMatrix(taxis,plusmin,corrmatrix)
end
%% Compare the contributions of all neurons to Jaw modes found at different time-points
for gg = 1:length(meta)
    lateDelWeights = mean(abs(blah(gg).tsModes(:,e1:e2)),2);        % Get the average weight (absolute value) for each neuron in the jaw modes found during the late delay period
    %RespWeights = mean(abs(blah(gg).tsModes(:,800:900)),2);        % Late response period
    MoveInitWeights = mean(abs(blah(gg).tsModes(:,511:601)),2);     % During movement initiation/very early response period
    
    %plotNeuronContributionScatter(lateDelWeights,MoveInitWeights)
end
%%
function plotTSCorr(tscorr,kinfeat,taxis)
figure();
imagesc(taxis(7:end-3),1:size(tscorr,1),tscorr);
colorbar
colormap('jet')
xlabel('Time since go-cue (s)')
ylabel('Cells')
plottitle = strcat('Time-varying correlation with',kinfeat);
title(plottitle)
end  % plotTSCorr

function plotCorrMatrix(taxis,plusmin,corrmatrix)
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
xlabel('Time since go-cue (s)','FontSize',15)
ylabel('Time since go-cue (s)','FontSize',15)
end

function plotNeuronContributionScatter(lateDelWeights,MoveInitWeights)
figure();
xx = linspace(0,max(max(lateDelWeights),max(MoveInitWeights)));
yy = xx;
scatter(lateDelWeights,MoveInitWeights,'filled')                % See if there is any relationship between contribution of neurons to delay per Jaw Mode vs Move initation Jaw Mode
hold on; plot(xx,yy,'LineStyle','--')
xlabel('Late Delay')
ylabel('Move initiation')
title('Average jaw mode weights')
end