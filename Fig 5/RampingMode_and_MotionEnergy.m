% Correspond ramping mode on single trials to motion energy
clear,clc,close all

% Add paths
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig3')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 6';
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
%meta = loadJEB13_ALMVideo(meta,datapth);
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
%% Find the Null and Potent Space and Calculate CDs

clearvars -except obj meta params me sav datapth kin

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -- Calculate coding directions from null and potent spaces --
% -----------------------------------------------------------------------
nSessions = numel(meta);
for sessix = 1:nSessions

    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat,obj(sessix));

    % -- null and potent spaces
    message = strcat('----Calculating N/P spaces for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    cond2use = [2 3 4 5]; % right hit, left hit, right miss, left miss
    cond2proj = [2:8,1];  % All conditions, last entry = all single trials
    nullalltime = 0; % use all time points to estimate null space if 1
    AWonly = 0;             % use only AW to find null and potent spaces 
    delayOnly = 0;          % use only delay period to find null and potent spaces
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj, nullalltime,AWonly,delayOnly);
    
    % -- coding dimensions
    message = strcat('----Calculating N/P CDs for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    cond2use = [1 2]; % right hits, left hits (corresponding to null/potent psths in rez)
    rampcond = 8;
    cond2proj = 1:6; % right hits, left hits, right miss, left miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
    cd_null(sessix) = getCDs_wRamping(rez(sessix).N_null_psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,rampcond);
    cd_potent(sessix) = getCDs_wRamping(rez(sessix).N_potent_psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,rampcond);
    
    cond2use = [2 3]; % right hits, left hits (corresponding to PARAMS.CONDITION)
    rampcond = 8;
    cond2proj = 2:5; % right hits, left hits, right miss, left miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
    cd_ramping(sessix) = getCDs_wRamping(obj(sessix).psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,rampcond);
end
%% Get projections onto Ramping mode on single trials
clearvars -except cd_null cd_potent cd_ramping meta obj params rez me

orthogonalize = 'non-orthog';                                       % Set to orthogonalize if you want the projections to be onto the orthogonalized CDs
disp('----Projecting single trials onto Ramping----')
cd = 'ramping';
cd_ramping = getSingleTrialProjs(cd_ramping,obj,cd,orthogonalize);
%%
% Specify time epoch by which you want to sort ME and ramping
delay = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent));
go = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent));
startix = find(obj(1).time>delay,1,'first');
stopix = find(obj(1).time<go,1,'last');

cond2use = 8;
ngroups = 4;
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

nTrials = length(sortix);
trixPerGroup = floor(nTrials/ngroups);                  % How many trials you want to be in each group
cnt = 1;
tempME = cell(1,ngroups);
tempRamp = cell(1,ngroups);
for g = 1:ngroups
    if g==ngroups
        ixrange = cnt:nTrials;
    else
        ixrange = cnt:(cnt+trixPerGroup);
    end
    tempME{g} = mean(sortedME(:,ixrange),2,'omitnan');
    tempRamp{g} = mean(sortedramp(:,ixrange),2,'omitnan');
    cnt = cnt+trixPerGroup+1;
end
grouped(sessix).ME = tempME;
grouped(sessix).ramp = tempRamp;


% figure();
% subplot(1,2,1)
% imagesc(obj(1).time,1:size(sortedME,2),sortedME')
% colorbar()
% %clim([0 100])
% % colormap('linspecer')
% title('Motion Energy')
% xlabel('Time from go cue (s)')
% ylabel('Trials')
% 
% subplot(1,2,2)
% imagesc(obj(1).time,1:size(sortedramp,2),sortedramp')
% colorbar()
% % colormap('linspecer')
% clim([5 40])
% title('Projection onto ramp mode')
% xlabel('Time from go cue (s)')
% 
% sgtitle([meta(sessix).anm ' ' meta(sessix).date])
% disp('hi')
end
%%
cols = {[0 0 0],[0.4 0.4 0.4],[0.6 0.6 0.6],[0.75 0.75 0.75]};
sm = 61;
for sessix = 1:length(meta)
    figure();
    for group = 1:ngroups
        subplot(2,1,1)
        toplot = mySmooth(grouped(sessix).ME{group},sm);
        plot(obj(sessix).time,toplot,'Color',cols{group},'LineWidth',1.5); hold on;
        title('Motion energy, grouped')
        ax1=gca;

        subplot(2,1,2)
        toplot = mySmooth(grouped(sessix).ramp{group},sm);
        plot(obj(sessix).time,toplot,'Color',cols{group},'LineWidth',1.5); hold on;
        title('Proj onto ramping mode')
        ax2 = gca;
    end
    legend(ax1,{'1st','2nd','3rd','4th'},'Location','best')
    xline(ax1,0,'k--','LineWidth',1)
    xlim(ax1,[-2.4 0.1])
    xlim(ax2,[-2.4 0.1])
    hold(ax1,'off');  hold(ax2,'off');
    sgtitle([meta(sessix).anm,' ',meta(sessix).date])
end
%%
cols = {[0 0 0],[0.4 0.4 0.4],[0.6 0.6 0.6],[0.75 0.75 0.75]};
sm = 61;
for sessix = 1:length(meta)
    figure();
    for group = 1:ngroups
        subplot(2,2,group)
        toplot = mySmooth(grouped(sessix).ME{group},sm);
        yyaxis left 
        plot(obj(sessix).time,toplot,'Color',cols{group},'LineWidth',1.5); hold on;
        toplot = mySmooth(grouped(sessix).ramp{group},sm);
        ylabel('Motion energy')
        yyaxis right 
        plot(obj(sessix).time,toplot,'Color',cols{group},'LineWidth',1.5,'LineStyle','--');
        title(['Group ',' ',num2str(group)])
        ylabel('Ramping mode')

        xline(0,'k--','LineWidth',1)
        xlim([-2.4 0.1])
        ax1=gca;
    end
    legend(ax1,{'Motion Energy','Ramping mode'},'Location','best') 
    sgtitle([meta(sessix).anm,' ',meta(sessix).date])
end
%%
cd_null_all = concatRezAcrossSessions(cd_null);
cd_potent_all = concatRezAcrossSessions(cd_potent);

%% Project single trials onto Null and Potent CDs
disp('----Projecting single trials onto CDlate----')
cd = 'late';

[cd_null,cd_potent] = getNPSingleTrialProjs(obj,cd,cd_null,cd_potent,rez);
%%
% -----------------------------------------------------------------------
% -- Coding Dimensions --
% -----------------------------------------------------------------------

spacename = 'Null';
figure();
%plotCDProj_Ramping(cd_null_all,cd_null,spacename)

spacename = 'Potent';
figure();
%plotCDProj_Ramping(cd_potent_all,cd_potent,spacename)

plotNP_CD_Ramp(cd_null_all,cd_null,cd_potent_all,cd_potent)