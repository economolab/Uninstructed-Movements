% From Inagaki 2019 Nature (Discrete Attractor Dynamics paper)
% "To compute selectivity, we took the spike rate difference between two correct trial types for each selective neuron.
% For peri-stimulus time histograms (Figs. 5, 6, Extended Data Fig. 8k), only correct trials were included. For the peri-stimulus time histograms and selectivity of the 
% random delay task (Fig. 6, Extended Data Fig. 8k), only spikes before the go cue were pooled. Spikes were averaged over 100 ms with a 1-ms sliding window."

clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Hazarded Delay');
addpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\functions\HazDel funcs');

% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for                            % all trials
params.condition(1) = {'R&hit&~stim.enable&~autowater&~early'};             % R 2AFC hits, no stim
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % L 2AFC hits, no stim

params.tmin = -1.6;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 31;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% vid features to use
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)
params.advance_movement = 0;

% Haz delay params
params.delay(1) = 0.3000;
params.delay(2) = 0.6000;
params.delay(3) = 1.2000;
params.delay(4) = 1.8000;
params.delay(5) = 2.4000;
params.delay(6) = 3.6000;
%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

meta = [];

% --- ALM --- 
meta = loadJEB11_ALMVideo(meta,datapth);
meta = loadJEB12_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written
%% LOAD DATA
% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);
%% For each session--Trials separated by delay length
% Get the trialIDs corresponding to each delay length
% Find the PSTH for R and L trials of each delay length
for sesh = 1:length(meta)
    met = meta(sesh);

    anm = obj(sesh).pth.anm;                  % Animal name
    date = obj(sesh).pth.dt;                  % Session date
    probenum = string(met.probe);       % Which probe was used

    del(sesh).delaylen = obj(sesh).bp.ev.goCue - obj(sesh).bp.ev.delay;       % Find the delay length for all trials
    conditions = {1,2};
    del(sesh).del_trialid = getDelayTrix(params(sesh),conditions,del(sesh));     % Group the trials in each condition based on their delay length
    del(sesh).delPSTH = getPSTHbyDel(params(sesh),del(sesh),obj(sesh));          % Get avg PSTH for each delay length
end
%% Get the selectivity for every neuron during each delay period length
smooth = 50;
for sesh = 1:length(meta)
    temp = cell(1,length(params(1).delay));
    for d = 1:length(params(1).delay)
        blah = mySmooth(del(sesh).delPSTH.right{d},smooth)-mySmooth(del(sesh).delPSTH.left{d},smooth);
        temp{d} = sqrt(blah.^2);
        %temp{d} = mySmooth(blah,smooth);
    end
    del(sesh).selectivity = temp;
end
%% Pool all neurons across sessions
sel_all = cell(1, length(params(1).delay));
for d = 1:length(params(1).delay)
    temp = [];
    for sesh = 1:length(meta)
        temp = [temp,del(sesh).selectivity{d}];
    end
    sel_all{d} = temp;
end 
%%
% smooth = 51;
% sesh = 1; del2use = 3;
% ncells = size(del(sesh).delPSTH.right{del2use},2);
% for i = 1:ncells
%     subplot(1,2,1)
%     R = mySmooth(del(sesh).delPSTH.right{del2use}(:,i),smooth);
%     L = mySmooth(del(sesh).delPSTH.left{del2use}(:,i),smooth);
%     plot(taxis,R);
%     hold on; plot(taxis,L);
%     xlim([-1.56 1.2])
%     hold off;
%     subplot(1,2,2)
%     plot(taxis,R-L);
%     xlim([-1.56 1.2])
%     pause
% end
% blah = mean(sel_all{del2use}(186:239,:),1);
% histogram(blah,20)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------Static Delay-----------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PARAMETERS
ctrlparams.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
ctrlparams.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
ctrlparams.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

ctrlparams.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for                            % all trials
ctrlparams.condition(1) = {'R&hit&~stim.enable&~autowater&~early'};             % R 2AFC hits, no stim
ctrlparams.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % L 2AFC hits, no stim

ctrlparams.tmin = -1.6;
ctrlparams.tmax = 2.5;
ctrlparams.dt = 1/200;

% smooth with causal gaussian kernel
ctrlparams.smooth = 31;

% cluster qualities to use
ctrlparams.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% vid features to use
ctrlparams.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

ctrlparams.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
ctrlparams.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)
ctrlparams.advance_movement = 0;

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

ctrlmeta = [];

% --- ALM --- 
ctrlmeta = loadJEB6_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB7_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadEKH1_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadEKH3_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJGR2_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJGR3_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB14_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB15_ALMVideo(ctrlmeta,datapth);

ctrlparams.probe = {ctrlmeta.probe}; % put probe numbers into params, one entry for element in ctrlmeta, just so i don't have to change code i've already written
%% LOAD DATA
% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[ctrlobj,ctrlparams] = loadSessionData(ctrlmeta,ctrlparams);
%% Get the selectivity for every neuron during each delay period length
smooth = 50;
for sesh = 1:length(ctrlmeta)
    blah = ctrlobj(sesh).psth(:,:,1)-ctrlobj(sesh).psth(:,:,2);
    temp = mySmooth(sqrt(blah.^2),smooth);
    ctrl(sesh).selectivity = temp;
end
%% Pool all neurons across sessions
temp = [];
for sesh = 1:length(ctrlmeta)
    temp = [temp,ctrl(sesh).selectivity];
end 
ctrl_sel_all = temp;
%% Take and plot the average selectivity across all neurons
taxis = ctrlobj(1).time;
del2use = 3;
nNeurons= size(sel_all{del2use},2);
col = [0 0 0];
alpha = 0.2;
toplot = mean(sel_all{del2use},2,'omitnan');            % Take the mean squared-selectivity across all neurons (Haz del)
err = 1.96*(std(sel_all{del2use},0,2,'omitnan')/nNeurons);
ax = gca;
shadedErrorBar(taxis, toplot, err ,{'Color',col,'LineWidth',2.5}, alpha, ax); hold on;

nNeurons= size(ctrl_sel_all,2);
col = [0.45 0.45 0.45];
alpha = 0.2;
toplot = mean(ctrl_sel_all,2,'omitnan');               % Take the mean squared-selectivity across all neurons (static del)
err = 1.96*(std(ctrl_sel_all,0,2,'omitnan')/nNeurons);
ax = gca;
shadedErrorBar(taxis, toplot, err ,{'Color',col,'LineWidth',1.5,'LineStyle','-.'}, alpha, ax); hold on;
xlim([-1.5 1.2])
xline(0,'LineStyle','--','Color','black')
xline(-1.3, 'LineStyle','--','Color','black')
xlabel('Time from delay onset (s)')
ylabel('Selectivity (spikes per s)')
legend('Hazarded delay, 1.2s trials','Fixed delay, 0.9s','Location','best')
