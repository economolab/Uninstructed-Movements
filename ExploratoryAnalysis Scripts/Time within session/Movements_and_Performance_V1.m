% Relating behavioral performance to uninstructed movements
clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
addpath(genpath(fullfile(utilspth,'Context_funcs')));
otherpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements';
addpath(genpath(fullfile(otherpth,'Decoding Analysis')));
addpath(genpath(fullfile(otherpth,'Context Switching\funcs')));
addpath(genpath(fullfile(otherpth,'functions')));

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~no&~stim.enable&~autowater&~early'};  % R hit 2AFC
params.condition(end+1) = {'L&hit&~no&~stim.enable&~autowater&~early'};  % L hit 2AFC
params.condition(end+1) = {'R&miss&~no&~stim.enable&~autowater&~early'};  % R hit 2AFC
params.condition(end+1) = {'L&miss&~no&~stim.enable&~autowater&~early'};  % L hit 2AFC

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


% params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
%     {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;

params.bctype = 'reflect'; % options are : reflect  zeropad  none
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
meta = loadJEB14_ALMVideo(meta,datapth);
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

%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = ['----Getting kinematic data for session' ,num2str(sessix), ' out of ' ,num2str(nSessions),'----'];
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
% %% Sanity check
% nfeats = length(kin(1).featLeg);
% cond2use = [2,3];
% figure();
% alpha = 0.2;
% for sessix = 1:numel(meta)
%     for feat = 1:nfeats
%         for c = cond2use
%             if c==2
%                 col = 'blue';
%             else
%                 col = 'red';
%             end
%             condtrix = params(sessix).trialid{c};
%             ax = gca;
%             %temperr = 1.96*(std(kin(sessix).dat(:,condtrix,feat),0,2,'omitnan')/length(condtrix));
%             plot(obj(1).time, kin(sessix).dat(:,condtrix,feat),'Color',col)
%             %shadedErrorBar(obj(1).time, mean(kin(sessix).dat(:,condtrix,feat),2,'omitnan'), temperr ,{'Color',col,'LineWidth',2}, alpha, ax)
%             hold on;
%         end
%         xline(-2.2,'LineStyle','--')
%         title(kin(1).featLeg{feat})
%         hold off;
%         pause
%     end
% end
%% Calculate all CDs and find single trial projections
clearvars -except obj meta params me sav kin

feats2use = [3,4,7,8,11:24,27,28,31,32,35,36,39,40:52,55];
disp('----Calculating MoveCDlate Mode----')
cond2use = [2,3];
cond2proj = [2,3,4,5];
regr = getMovement_CDlate(obj,params,cond2use,cond2proj,kin,feats2use);

disp('----Projecting single trials of movement onto MoveCDlate----')
cd = 'MOVECDlate';
regr = getMove_SingleTrialProjs(regr,obj,cd,kin,feats2use);

%% For single sessions, show all trials projected onto MoveCDlate, and condition-averaged projections
colors = {'blue','red',[0.5 0.5 1],[1 0.5 0.5]};
alpha = 0.2;
for sessix = 1:length(meta)
    figure();
    for c = 1:length(cond2proj)
        condtrix = params(sessix).trialid{cond2proj(c)};
        temp = regr(sessix).singleMoveCDProj(:,condtrix);
        if ~isempty(temp)
            %plot(obj(1).time,temp,'Color',colors{c})
            regr(sessix).condavgProj(:,c) = mean(temp,2,'omitnan');
            temperr = 1.96*(std(temp,0,2,'omitnan')/length(condtrix));
            ax = gca;
            shadedErrorBar(obj(1).time, mean(temp,2,'omitnan'), temperr ,{'Color',colors{c},'LineWidth',2}, alpha, ax)
            hold on;
            xline(-2.2,'Linestyle','--','Color','black')
        else
            regr(sessix).condavgProj(:,c) = NaN(length(obj(1).time),1);
        end
    end
end
%% Get condition-averaged MoveCDContext across all sessions
% Take avg across sessions 
sessavg = NaN(length(obj(1).time),length(cond2proj));        % (time x conditions)
sessstd = NaN(length(obj(1).time),length(cond2proj));
for c = 1:length(cond2proj)
    temp = [];
    for sessix = 1:length(meta)
        temp = [temp,regr(sessix).condavgProj(:,c)];
    end
    sessavg(:,c) = mean(temp,2,'omitnan');
    sessstd(:,c) = std(temp,0,2,'omitnan');
end
%% Plot condition-averaged MoveCDlate across all sessions
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
delay = median(obj(1).bp.ev.delay)-median(obj(1).bp.ev.(params(1).alignEvent));
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
% Set plot params and plot
colors = {'blue','red',[0.5 0.5 1],[1 0.5 0.5]};
alpha = 0.2;
figure();
cond2proj = [1,2,3,4];
for c = 1:length(cond2proj)
    temperr = 1.96*(sessstd(:,c)/length(meta));
    toplot = sessavg(:,c);
    toplot(isinf(toplot)) = NaN;
    toplot = fillmissing(toplot,'previous');
    ax = gca;
    shadedErrorBar(obj(1).time, toplot, temperr ,{'Color',colors{c},'LineWidth',2}, alpha, ax); hold on;
end
xline(0,'k','Linestyle','--')
xline(delay,'k','Linestyle','-.')
xline(samp,'k','Linestyle','-.')
xlabel('Time from go cue (s)')
ylabel('Projection onto Movement-CDLate (a.u.)')
xlim([trialstart 2.5])
%%
figure();
start = find(obj(1).time>-0.9,1,'first');
stop = find (obj(1).time<-0.05,1,'last');
for sessix = 1:length(meta)
    temp = regr(sessix).singleMoveCDProj;
    temp = mean(temp(start:stop,:),1,'omitnan');
    plot(temp)
    pause
end
%% Sort Movement-CDlate according to value in late delay period
% Remove the first and last 'X' trials from analysis
nBufferTrix = 20;
for sessix = 1:length(meta)
    temp = regr(sessix).singleMoveCDProj(:,nBufferTrix:end-nBufferTrix);
end

%% FUNCTIONS
