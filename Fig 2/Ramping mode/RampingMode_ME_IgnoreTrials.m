% Correspond ramping mode on single trials to motion energy
clear,clc,close all

% Add paths
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig3')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 5';
addpath(genpath(fullfile(figpth,'funcs')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 2';
addpath(genpath(fullfile(figpth,'funcs')));
%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off
params.condition(end+1) = {'R&no&~stim.enable&~autowater&~early'};              % no right, no stim, aw off
params.condition(end+1) = {'L&no&~stim.enable&~autowater&~early'};              % no left, no stim, aw off
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % all hits, no stim, aw off

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
% params.quality = {'Excellent','Great','Good','Fair','Multi'};

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;
%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

meta = [];

% --- ALM --- 
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB13_ALMVideo(meta,datapth);
%meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written
%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);

% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end
%%
clearvars -except obj meta params me sav datapth kin

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials

nSessions = numel(meta);
for sessix = 1:nSessions
    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat,obj(sessix));
    cond2use = [2 3 6 7]; % right hits, left hits (corresponding to PARAMS.CONDITION)
    rampcond = 8;
    cond2proj = 2:7;  % right hits, left hits, right miss, left miss, right no, left no (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
    cd_ramping(sessix) = getCDs_wRamping(obj(sessix).psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,rampcond);
end
%% Get projections onto Ramping mode on single trials
clearvars -except cd_ramping meta obj params me

orthogonalize = 'non-orthog';                                       % Set to orthogonalize if you want the projections to be onto the orthogonalized CDs
disp('----Projecting single trials onto Ramping----')
cd = 'ramping';
cd_ramping = getSingleTrialProjs(cd_ramping,obj,cd,orthogonalize);
%%
colors = getColors_Updated();
cond2proj = [1 2 5 6];
mecond2proj = [2 3 6 7];
for sessix = 1:length(meta)
    for c = 1:length(cond2proj)
        switch c
            case 1
                col = colors.rhit;
                linestyl = '--';
                lw = 1.5;
            case 2
                col = colors.lhit;
                linestyl = '--';
                lw = 1.5;
            case 3
                col = colors.rhit;
                linestyl = '-';
                lw = 2;
            case 4
                col = colors.lhit;
                linestyl = '-';
                lw = 2;
        end
        ax1 = subplot(1,2,1);
        toplot = cd_ramping(sessix).cd_proj(:,cond2proj(c),4);
        plot(obj(sessix).time,mySmooth(toplot,31),'Color',col,'LineStyle',linestyl,'LineWidth',lw)
        hold(ax1,'on');

        ax2 = subplot(1,2,2);
        condME = me(sessix).data(:,params(sessix).trialid{mecond2proj(c)});
        condME = mean(condME,2,'omitnan');
        plot(obj(sessix).time,mySmooth(condME,31),'Color',col,'LineStyle',linestyl,'LineWidth',lw)
        hold(ax2, 'on');
    end
    legend(ax1,'R hit', 'L hit','R ignore','L ignore')
    title(ax1,'Ramping')
    title(ax2,'ME')
    hold (ax1, 'off'); hold (ax2, 'off');
    pause
end

%%
% Specify time epoch by which you want to sort ME and ramping
delay = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent));
go = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent));
startix = find(obj(1).time>delay,1,'first');
stopix = find(obj(1).time<go,1,'last');

cond2use = 8;
for sessix = 1:length(meta)
% Get trials for 2AFC
trix2use = params(sessix).trialid{cond2use};

% Get ME and projections onto Ramping Mode for these trials
MEtrix = me(sessix).data(:,trix2use);
rampingtrix = cd_ramping(sessix).singleProj(:,trix2use);

% Find avg ME during delay on each trial
avgME = mean(MEtrix(startix:stopix,:),1,'omitnan');
% Sort ME and Ramping proj by ME descending order
[~,sortix] = sort(avgME,'descend');
sortedME = MEtrix(:,sortix);
sortedramp = rampingtrix(:,sortix);

delayME = [];
delayramp = [];
for trix = 1:size(MEtrix,2)
    delayME = [delayME,MEtrix(startix:stopix,trix)];
    delayramp = [delayramp,rampingtrix(startix:stopix,trix)];
end
R2 = corrcoef(delayME,delayramp);
R2 = R2(2);
obj(sessix).R2 = R2;

if sessix == 1
figure();
subplot(1,2,1)
imagesc(obj(1).time,1:size(sortedME,2),sortedME')
colorbar()
clim([0 80])
colormap('linspecer')
title('Motion Energy')
xlabel('Time from go cue (s)')
ylabel('Trials')
xlim([-2.5 0])

subplot(1,2,2)
imagesc(obj(1).time,1:size(sortedramp,2),sortedramp')
colorbar()
colormap('linspecer')
clim([0 40])
title('Projection onto ramp mode')
xlim([-2.5 0])
xlabel('Time from go cue (s)')
sgtitle([meta(sessix).anm ' ' meta(sessix).date '  ; R2 = ' num2str(R2)])
end
end
