% DECODING CDlate FROM ALL KINEMATIC FEATURES
clear,clc,close all

whichcomp = 'LabPC';                                                % LabPC or Laptop

% Base path for code depending on laptop or lab PC
if strcmp(whichcomp,'LabPC')
    basepth = 'C:\Code';
elseif strcmp(whichcomp,'Laptop')
    basepth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code';
end

% add paths
utilspth = [basepth '\Munib Uninstruct Move\uninstructedMovements_v2'];
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
figpth = [basepth  '\Uninstructed-Movements\Fig 2'];
addpath(genpath(fullfile(figpth,'funcs')));
%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};             % left hits, no stim, aw off

params.tmin = -3;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;
params.bctype = 'reflect'; % options are : reflect  zeropad  none
%% SPECIFY DATA TO LOAD

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written
%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);
%%
% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    disp(['Loading ME for session ' num2str(sessix)])
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% Example session of tongue trajectories for bottom cam (across contexts)
sess2use = 1;
cond2use = 2:5;
feat2use = {'left_tongue'};
featix = find(strcmp(obj(1).traj{1}(1).featNames,feat2use));

go = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent));
goix = find(obj(1).time<go,1,'last');
resp = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent))+2.5;
respix = find(obj(1).time<resp,1,'last');

colors = getColors();
nTrixPlot = 14;
for sessix = sess2use
    figure();
    for feat = 1:length(feat2use)
        for c = 1:length(cond2use)
            switch c
                case 1
                    col = colors.rhit;
                case 2
                    col = colors.lhit;
                case 3
                    col = colors.rhit_aw;
                case 4
                    col = colors.lhit_aw;
            end
            cond = cond2use(c);
            condtrix = params(sessix).trialid{cond};
            randtrix = randsample(condtrix,nTrixPlot);
            for trix = 1:nTrixPlot
                currtrix = randtrix(trix);
                condfeat = squeeze(obj(sessix).traj{1}(currtrix).ts(:,1:2,featix));
                plot(condfeat(:,1),condfeat(:,2),'Color',col); hold on
            end
        end
    end
    title(['Tongue traj, side view; ' meta(sessix).anm meta(sessix).date])
end
%% Example session of tongue trajectories for bottom cam (across contexts)
sess2use = 2;
cond2use = 2:5;
feat2use = {'top_tongue'};
featix = find(strcmp(obj(1).traj{2}(1).featNames,feat2use));

go = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent));
goix = find(obj(1).time<go,1,'last');
resp = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent))+2.5;
respix = find(obj(1).time<resp,1,'last');

colors = getColors();
nTrixPlot = 14;
for sessix = sess2use
    figure();
    for feat = 1:length(feat2use)
        for c = 1:length(cond2use)
            switch c
                case 1
                    col = colors.rhit;
                case 2
                    col = colors.lhit;
                case 3
                    col = colors.rhit_aw;
                case 4
                    col = colors.lhit_aw;
            end
            cond = cond2use(c);
            condtrix = params(sessix).trialid{cond};
            randtrix = randsample(condtrix,nTrixPlot);
            for trix = 1:nTrixPlot
                currtrix = randtrix(trix);
                condfeat = squeeze(obj(sessix).traj{2}(currtrix).ts(:,1:2,featix));
                plot(-mySmooth(condfeat(:,1), 5),-mySmooth(condfeat(:,2), 5),'Color',col); hold on
            end
        end
    end
    title(['Tongue traj, bottom view; ' meta(sessix).anm meta(sessix).date])
end
%% Get tongue length and angle for individual licks for all animals and sessions
% Stored in variable 'tonguefeat': (1 x nSessions)
% Each session will have fields 'tongue_length' and 'tongue_angle'
% Each field will have subfields 'RAFC', 'LAFC','RAW','LAW'

nLicks = 5;
cond2use = 2:5;
condfns = {'RAFC','LAFC','RAW','LAW'};
feat2use = {'tongue_angle','tongue_length'};
featix = NaN(1,length(feat2use));
for f = 1:length(feat2use)
    currfeat = feat2use{f};
    currix = find(strcmp(kin(1).featLeg,currfeat));
    featix(f) = currix;
end

% Get times for which you want to look at the tongue
go = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent));
goix = find(obj(1).time<go,1,'last');
resp = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent))+2.5;
respix = find(obj(1).time<resp,1,'last');

for sessix = 1:length(meta)                                                 % For every session...
    for feat = 1:length(feat2use)                                           % For tongue length and angle...
        metric = feat2use{feat};
        for c = 1:length(cond2use)                                          % For each condition...
            cond = cond2use(c);
            condtrix = params(sessix).trialid{cond};                        % Get the trials for this condition
            condfeat = squeeze(kin(sessix).dat(:,condtrix,featix(feat)));   % Get the kinematic data for these trials and current feature
            respfeat = condfeat(goix:respix,:);                             % Take only desired post-goCue time points

            nTrials = size(respfeat,2);                                     % number of trials in this condition
            nTimePts = size(respfeat,1);                                    % number of time-points

            lickDur = 10;
            temp = NaN(lickDur,nLicks,nTrials);                                    % Pre-allocate (nLicks x nTrials) for this condition
            avgDur = NaN(nLicks, nTrials);
            for trix = 1:nTrials                                            % For each trial in this condition...
                licks.start = [];
                licks.stop = [];
                % Find indices corresponding to lick start and lick end
                for pt = 4:nTimePts-4                                       % For every time point...
                    % If the tongue is not visible at the current time
                    % point but is visible for the following three time
                    % pts...
                    if isnan(respfeat(pt,trix))&&~isnan(respfeat(pt+1,trix))&&~isnan(respfeat(pt+2,trix))...
                            &&~isnan(respfeat(pt+3,trix))
                        % Take the subsequent time-point to be the start of a lick
                        licks.start = [licks.start,pt+1];
                        % If the tongue is visible at the current and previous two time
                        % points but is not visible for the next time point...
                    elseif ~isnan(respfeat(pt,trix))&&~isnan(respfeat(pt-1,trix))&&...
                            ~isnan(respfeat(pt-2,trix))&&isnan(respfeat(pt+1,trix))
                        % Take this time point to be the end of a lick
                        licks.stop = [licks.stop,pt];
                    end
                end


                for lix = 1:nLicks                                                  % For each lick that you want to look at
                    if lix<size(licks.start,2)&&lix<size(licks.stop,2)              % If there are this many licks in the current trial...
                        currstart = licks.start(lix);
                        currstop = licks.stop(lix);
                        if currstart<currstop
                            currtongue = respfeat(currstart:currstop,trix);             % Take the tongue length/angle for this current lick and store it
                            currdur = size(currtongue,1);                               % Get the duration in frames of each individual lick
                            avgDur(lix,trix) = currdur;
                            tongueinterp = interp1(1:length(currtongue),currtongue, linspace(1,length(currtongue),lickDur));
                            temp(:,lix,trix) = tongueinterp';
                        end
                    end
                end
            end
            tonguefeat(sessix).(metric).(condfns{c}) = temp;
            tonguefeat(sessix).lickDur.(condfns{c}) = avgDur;
        end
    end
end
clearvars -except tonguefeat obj params meta kin me cond2use condfns feat2use nLicks
%% Find average tongue length and angle for each condition within each session
feat2use = {'tongue_length','tongue_angle','lickDur'};
for feat = 1:length(feat2use)
    for c = 1:length(cond2use)
        tempAll = [];
        varAll = [];
        for sessix = 1:length(meta)
            currtongue = tonguefeat(sessix).(feat2use{feat}).(condfns{c});
            if ~strcmp(feat2use{feat},'lickDur')
                tempmean = squeeze(mean(currtongue,3,'omitnan'));
                tempvar = var(currtongue,0,3,'omitnan');
                tempvar = mean(tempvar,1,'omitnan');
                tempAll = cat(3,tempAll,tempmean);
                varAll = [varAll,tempvar'];
            else
                tempmean=squeeze(mean(currtongue,2,'omitnan'));
                tempAll = [tempAll, tempmean];
            end
        end
        avgtongue.(feat2use{feat}).(condfns{c}) = tempAll;
        vartongue.(feat2use{feat}).(condfns{c}) = varAll;
    end
end
%% Plot average tongue angle for 1st, 2nd, 3rd, etc. licks for each session
feat2use = {'tongue_angle'};
condfns = {'RAFC','RAW','LAFC','LAW'};
lickDur = 10;
alph = 0.2;
colors = getColors;
for sessix = 1:length(meta)
    figure();
for feat = 1:length(feat2use)
    for lix = 1:nLicks
        subplot(1,nLicks,lix)
        for c = 1:length(condfns)
            switch c
                case 1
                    col = colors.rhit;
                case 2
                    col = colors.rhit_aw;
                case 3
                    col = colors.lhit;
                case 4
                    col = colors.lhit_aw;
            end
            currtongue = tonguefeat(sessix).(feat2use{feat}).(condfns{c})(:,lix,:);
            nTrials = size(currtongue,3);
            toplot = squeeze(mean(currtongue,3,'omitnan'));
            err = 1.96*(std(currtongue,0,3,'omitnan')./sqrt(nTrials));
            ax = gca;
            shadedErrorBar(1:lickDur,toplot,err,{'Color',col,'LineWidth',2},alph, ax); hold on;
        end
        title(['Lick ' num2str(lix)])
        if lix==1&&strcmp(feat2use{feat},'tongue_length')
            ylabel('Length')
        elseif lix==1&&strcmp(feat2use{feat},'tongue_angle')
            ylabel('angle')
        end
        xlim([0 11])
    end
    sgtitle([feat2use{feat}])
end
end
%% Plot average tongue length for 1st, 2nd, 3rd, etc. licks across sessions
feat2use = {'tongue_length'};
condfns = {'RAFC','RAW','LAFC','LAW'};
lickDur = 10;
alph = 0.2;
colors = getColors;
for sessix = 1:length(meta)
    figure();
    for feat = 1:length(feat2use)
        for lix = 1:nLicks
            for c = 1:length(condfns)
                switch c
                    case 1
                        col = colors.rhit;
                        sub = lix;
                    case 2
                        col = colors.rhit_aw;
                        sub = lix;
                    case 3
                        col = colors.lhit;
                        sub = lix+nLicks;
                    case 4
                        col = colors.lhit_aw;
                        sub = lix+nLicks;
                end
                subplot(2,nLicks,sub)
                currtongue = tonguefeat(sessix).(feat2use{feat}).(condfns{c})(:,lix,:);
                nTrials = size(currtongue,3);
                toplot = squeeze(mean(currtongue,3,'omitnan'));
                err = 1.96*(std(currtongue,0,3,'omitnan')./sqrt(nTrials));
                ax = gca;
                shadedErrorBar(1:lickDur,toplot,err,{'Color',col,'LineWidth',2},alph, ax); hold on;
                xlim([0 11])
            end
            title(['Lick ' num2str(lix)])
            if lix==1&&strcmp(feat2use{feat},'tongue_length')
                ylabel('Length')
            elseif lix==1&&strcmp(feat2use{feat},'tongue_angle')
                ylabel('angle')
            end
        end
        sgtitle([feat2use{feat}])
    end
end
%% Plot average duration of each lick (1st, 2nd, 3rd) across sessions and animals
camrate = 400;
framedur = 1/camrate;

sigcutoff = 0.05;
feat2use = {'lickDur'};
cond2use = {'RAFC','RAW','LAFC','LAW'};
for feat = 1:length(feat2use)
    figure();
    for lix = 1:nLicks
        temp = [];
        for c = 1:length(cond2use)
            cond = cond2use{c};
            curr = avgtongue.(feat2use{feat}).(cond)(lix,:)';
            temp = [temp, curr];
        end
        temp = framedur*temp*1000;           % Go from lick duration in frames to ms

        subplot(1,nLicks, lix)
        for c = 1:length(cond2use)
            switch c
                case 1
                    col = colors.rhit;
                case 2
                    col = colors.rhit_aw;
                case 3
                    col = colors.lhit;
                case 4
                    col = colors.lhit_aw;
            end
            bar(c,mean(temp(:,c),1,'omitnan'),'FaceColor',col,'EdgeColor',col); hold on
            xx = c*ones(length(meta),1);
            scatter(xx,temp(:,c),15,'filled','MarkerFaceColor','black')

            for test = 1:(length(cond2use)/2)
                x = temp(:,(test*2)-1);
                y = temp(:,test*2);
                hyp = ttest(x,y,'Alpha',sigcutoff);
                if hyp&&test==1
                    scatter(1.5,20,20,'*','MarkerEdgeColor','black')
                elseif hyp&&test==2
                    scatter(3.5,20,20,'*','MarkerEdgeColor','black')
                end
            end
            xlim([0 5])
            ylim([0 21])
            if lix==1
                ylabel('Duration (ms)')
            end
        end
        title(['Lick ' num2str(lix)])
    end
    sgtitle([feat2use{feat}])
end