clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')));
rmpath(genpath(fullfile(utilspth,'mc_stim/')));

% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM --- 
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

% meta = meta(1:2);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

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
disp('Loading Motion Energy')
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end


%% 

clearvars -except obj meta params me sav

dt = params(1).dt;

maxlag = int64(5 ./ dt);

for sessix = 1:numel(meta)

    % -- neural data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat);

    % -- motion energy
    me_zscored = zscore(me(sessix).data); % me(sessix).data ./ max(me(sessix).data);

    % for each neuron
    for cluix = 1:size(trialdat_zscored,3)
        % for each trial
        for trix = 1:size(trialdat_zscored,2)
            % cross correlation
            N = trialdat_zscored(:,trix,cluix);
%             [r{sessix}(:,trix,cluix),lags] = xcorr(N,me_zscored(:,trix), maxlag, 'normalized');
            [r{sessix}(:,trix,cluix),lags] = xcorr(me_zscored(:,trix),N, maxlag, 'normalized');
        end

    end

    
end

r_mean_across_trials = cellfun(@(x) squeeze(nanmean(x,2)),r,'UniformOutput',false);
r_mean_across_neurons = cellfun(@(x) squeeze(nanmean(x,2)),r_mean_across_trials,'UniformOutput',false);
rr = cell2mat(r_mean_across_neurons);

%% plot lags against neurons for trial-averaged data

t = lags .* dt;
ix1 = find(t>=-0.1,1,'first');
ix2 = find(t<=0.1,1,'last');

f = figure;
t = tiledlayout('flow');
for i = 1:numel(r_mean_across_trials)
    ax = nexttile;
    temp = r_mean_across_trials{i};
    [~,ix] = sort(mean(temp(ix1:ix2,:),1));
    
    imagesc(lags .*dt ,1:size(temp,2),temp(:,ix)')
    title([meta(i).anm '_' meta(i).date],'Interpreter','none')
end
xlabel(t,'lags (s)')
ylabel(t,'neurons')
title(t,'cross correlation b/w motion energy and single trial fr')

%% mean lag for all neurons

close all


t = lags .* dt;

rr_mean = mean(rr,2,'omitnan');
rr_err = std(rr,[],2,'omitnan') ./ numel(meta);

[maxs,imaxs] = max(abs(rr)); 
tmaxs = t(imaxs);

meanalph = 0.5;
linealph = 0.8;
lw = 2;
cols = linspecer(size(rr,2));

f = figure; 
ax = gca;
% ax = subplot(2,1,1);
hold on;
for i = 1:size(rr,2)
    patchline(t,rr(:,i),'EdgeColor',cols(i,:),'LineWidth',1,'EdgeAlpha',linealph)
    scatter(tmaxs(i),rr(imaxs(i),i),60,'MarkerEdgeColor','none','MarkerFaceColor',cols(i,:))
end
shadedErrorBar(t,rr_mean,rr_err,{'Color','k','LineWidth',lw},meanalph,ax)
% xlim([-0.1 0.1])
xlabel('Time shift (s)')
ylabel('correlation')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ax = subplot(2,1,2);
% tt = tmaxs(tmaxs > -5);
% histogram(tmaxs,5,'EdgeColor','none','FaceColor','k','FaceAlpha',0.7);
% % xlim([-2.5 0.05])
% xlabel('Time shift (s)')
% ylabel('occurrences')








