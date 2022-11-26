clear,clc,close all

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\cd_code_jackie'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\kinematicModes'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Utils'));
rmpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context\Context_funcs'));
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


params.tmin = -3;
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
%%  Plot kinematic measurements
params.kinfind = 'vel';
for sessix = 1:length(meta)         % For all loaded sessions...
    clear orthproj orthModes proj kinmodes kinfns numTrials kin orthmode mode varexp
    
    obj = objs{sessix};
    met = meta(sessix);
    [anm,date,probenum,taxis] = getSessMeta(obj,met);
    
    %%% GET FEATURE KINEMATICS (for all trial types) %%%
    conditions = {1,2};     
    kin = struct();
    kin.MEinterp = getME(obj,met,params,taxis,conditions);                  % Find interpolated motion energy
    kinfeat = findAllFeatKin(params, taxis, obj,conditions,met);            % Find all kinematic features (jaw, nose, tongue vel)
    kin.jawVel = kinfeat.jawVel; kin.noseVel = kinfeat.noseVel; kin.tongueVel = kinfeat.tongueVel;
    %kin.tongueAngle = findTongueAngle(taxis, obj, met,params,conditions);   % Find tongue angle
   
    kinfns = fieldnames(kin);

    % Plot Kinematic Feature for every trial in chronological order
    featName = 'MEinterp';
    trialnums = [met.trialid{1}; met.trialid{2}];                                   % Get the trial numbers for all trials that fit the desired conditions
    triallabel = [zeros(length(met.trialid{1}),1);ones(length(met.trialid{2}),1)];      % Label all trials. 2AFC trials = 0; AW trials = 1
    [~,trialsortix] = sort(trialnums,'ascend');                                     % Sort all included trials into chronological order
    kinsorted(sessix).(featName) = kin.(featName)(:,trialsortix); kinsorted(sessix).trialsortix = trialsortix;    % Sort the kinematic features and trial labels into chronological order
    kinsorted(sessix).triallabel = triallabel(trialsortix);

    %%% Plot example trials from 2AFC and AW %%%
    conditions = {1,2};
    featName = 'MEinterp';
    kin_by_cond(sessix).(featName) = NaN(length(taxis),numel(conditions));
    cnt = 0;
    for c = 1:numel(conditions)
        le = length(met.trialid{conditions{c}});
        if c==1
            trixrange = 1:le;
        else
            trixrange = cnt+1:cnt+le;
        end
        kin_by_cond(sessix).(featName)(:,c) = mySmooth(mean(kin.(featName)(:,trixrange),2,'omitnan'),21);
        cnt = cnt+le;
    end

%     % Sort the trials by jaw velocity during the late delay period
%     % Find the average jaw velocity during specified time points (on each trial)
%     startix = find(taxis>=-0.7, 1, 'first');
%     stopix = find(taxis<=-0.05, 1, 'last');
%     val.right = nanmean(kin_by_cond.(featName){1}(startix:stopix, :), 1);
%     val.left = nanmean(kin_by_cond.(featName){2}(startix:stopix, :), 1);
%     nanix = find(isnan(val.right)); nanix = find(isnan(val.left));                   % Get rid of trials where jaw velocity is always NaN
%     val.right(nanix) = [];  val.left(nanix) = [];
% 
%     % Sort the average jaw velocities in descending order and save the trial
%     % order
%     sort_by_cond = cell(1,numel(conditions));
%     [sort_by_cond{1}, six1] = sort(val.right, 'descend'); kin_by_cond.(featName){1} = kin_by_cond.(featName){1}(:,six1);
%     [sort_by_cond{2}, six2] = sort(val.left, 'descend'); kin_by_cond.(featName){2} = kin_by_cond.(featName){2}(:,six2);
% 
%     figure();
%     numTrixPlot = 30;
%     rangetoPlot = 1:numTrixPlot;
%     tFirstLick = find(taxis==0);
%     for i=1:numel(conditions)
%         if i==1
%             subplot(1,2,1);
%         elseif i==2
%             subplot(1,2,2);
%         end
%         imagesc(taxis,1:numTrixPlot,kin_by_cond.(featName){i}(:,1:numTrixPlot)'); colormap("hot");
%         xlim([-2.6, 0])
%         xlabel('Time before firstLick (s)','FontSize',13)
%         ylabel('Trials','FontSize',13)
%         if i==1
%             title('2AFC trials')
%         elseif i==2
%             title('AW trials')
%         end
%         c=colorbar;
%         clim([0 0.5])
%         ylabel(c,featName,'FontSize',12,'Rotation',90);
%     end
%     figtitle =  strcat('Example Trials from',anm,date,' ;  ','Probe ',probenum);  % Name/title for session
%     sgtitle(figtitle,'FontSize',16)
end
%% Find Context Mode and project onto it %%%
for sessix = 1:length(meta)
    met = meta(sessix);
    % Find context mode
    ev = objs{sessix}.bp.ev;
    [rez(sessix).Context_mode] = findCDContext(objs,sessix,ev,params);
    cond2proj = [1,2];

    % Project onto mode for each context
    temp = NaN(length(taxis),2);
    for c = 1:length(cond2proj)
        temp(:,c) = mySmooth(objs{sessix}.psth(:,:,c)*rez(sessix).Context_mode,51);
    end
    rez(sessix).Context_latent = temp;                  % Save context latents as (time x condition)
    clear temp

    % Project each trial onto context mode
    condlatents = cell(1,length(cond2proj));
    for c = 1:length(cond2proj)
        triallatents = NaN(length(taxis),length(met.trialid{c}));
        for trix = 1:length(met.trialid{c})
            triallatents(:,trix) = mySmooth(objs{sessix}.trialpsth_cond{c}(:,:,trix)*rez(sessix).Context_mode,51);
        end
        condlatents{c} = triallatents;
    end
    temp = [condlatents{1},condlatents{2}];
    kinsorted(sessix).Context_latent = temp(:,kinsorted(sessix).trialsortix); 
end
%%
for sessix = 5%1:length(meta)
met = meta(sessix);
trialnums = [met.trialid{1}; met.trialid{2}];

figure();
subplot(1,2,1)
imagesc(taxis,1:size(kinsorted(sessix).(featName),2),kinsorted(sessix).(featName)')
colormap('parula')
tmax = -0.5;
xlim([-2.6 tmax])
c = colorbar; clim([0 6])
hold on;

Nplotted = 0;
triallabel = kinsorted(sessix).triallabel;
for j = 1:length(trialnums)
    Nplotted = Nplotted+1;
    if triallabel(j)==0                                                         % If a 2AFC trial, plot a black rectangle on the side
        plot([0 0]+(tmax-0.05), -0.5+Nplotted+[0.1 0.9], 'Color',[0.75 0.75 0.75], 'LineWidth', 10);
    end
    if triallabel(j)==1                                                         % If an AW trial, plot a magenta rectangle on the side
        plot([0 0]+(tmax-0.05), -0.5+Nplotted+[0.1 0.9], 'magenta', 'LineWidth', 10);
    end
end
ylabel(c,featName,'FontSize',12,'Rotation',90);

s = subplot(1,2,2);
imagesc(taxis,1:size(kinsorted(sessix).Context_latent,2),kinsorted(sessix).Context_latent')
tmax = -0.5;
xlim([-2.6 tmax])
colormap(s,flipud(parula))
c = colorbar;
clim([-15 5])
hold on;
Nplotted = 0;
for j = 1:length(trialnums)
    Nplotted = Nplotted+1;
    if triallabel(j)==0                                                         % If a 2AFC trial, plot a black rectangle on the side
        plot([0 0]+(tmax-0.05), -0.5+Nplotted+[0.1 0.9], 'Color',[0.75 0.75 0.75], 'LineWidth', 10);
    end
    if triallabel(j)==1                                                         % If an AW trial, plot a magenta rectangle on the side
        plot([0 0]+(tmax-0.05), -0.5+Nplotted+[0.1 0.9], 'magenta', 'LineWidth', 10);
    end
end
ylabel(c,'ContextMode','FontSize',12,'Rotation',90);
figtitle =  strcat('Full session',met.anm,met.date,' ;  ','Probe ',probenum);  % Name/title for session
sgtitle(figtitle,'FontSize',16)
ylabel('Trials')
xlabel('Time before firstLick (s)')
end
%%
% cond2plot = [1,2];
% colors = {[0 0 0],'magenta'};
% for sessix = 1:length(meta)
%     figtitle =  strcat('Avg ME and CDContext for',anm,date);  % Name/title for session
%     figure();
%     for c = 1:length(cond2plot)
%     subplot(2,1,1)
%     plot(taxis,kin_by_cond(sessix).MEinterp(:,c),'Color',colors{c},'LineWidth',2);
%     hold on;
%     ylabel('Motion Energy')
%     xlim([-2.9 -0.5])
% 
%     subplot(2,1,2)
%     plot(taxis,rez(sessix).Context_latent(:,c),'Color',colors{c},'LineWidth',2);
%     hold on;
%     xlabel('Time before firstLick (s)')
%     ylabel('Context Mode projection (a.u.)')
%     xlim([-2.9 -0.5])
%     sgtitle(figtitle,'FontSize',16)
%     end
% 
% end
%% FUNCTIONS

function MEinterp = getME(obj,met,params,taxis,conditions)
[met,mov,me] = assignEarlyTrials(obj,met,params);
[temp,~] = findInterpME(taxis,conditions, met,mov,me,params,obj);
MEinterp = [];
for c = 1:numel(conditions)
    MEinterp = [MEinterp, temp{c}];
end
end

%% PLOTTING FUNCTIONS
function plotModeFeat_SingleTrix(met,kinfns,taxis,anm,date,probenum,kin,params)
numTrials = length(met.trialid{1})+length(met.trialid{2});
l1 = length(met.trialid{1});
l2 = l1+length(met.trialid{2});

figure();
for i=1:numel(kinfns)
    % Kinematic feature tracking
    subplot(1,numel(kinfns),(i));
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
    xlab = strcat('Time since ',{' '},params.alignEvent,{' '},'(s)');
    xlabel(xlab)
    title(kinfns{i});
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