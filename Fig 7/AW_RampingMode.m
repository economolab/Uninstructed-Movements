% Find context selectivity across only Null and potent dimensions which are selective for context
% Take the difference between contexts for each selective dimension.  Take the mean
% across all selective dimensions within a session
% -------------------------------------------------------------------------------------
% Using all 2AFC and all AW trials to find the Null and Potent Spaces
% -------------------------------------------------------------------------------------
clear,clc,close all

% add paths
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig3')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 4';
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context switching')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 7';
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context_funcs')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 6';
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
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};                   % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};                   % error left, no stim, aw off
params.condition(end+1) = {'~early&~stim.enable&~autowater'};                          %  no stim, aw off

params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};             % right hits, no stim, aw on
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};             % left hits, no stim, aw on
params.condition(end+1) = {'R&miss&~stim.enable&autowater'};                   % error right, no stim, aw on
params.condition(end+1) = {'L&miss&~stim.enable&autowater'};                   % error left, no stim, aw on
params.condition(end+1) = {'~early&~stim.enable&autowater'};                          %  no stim, aw on

params.condition(end+1) = {'hit&~early&~stim.enable&~autowater'};                          %  hit, no stim, aw off
params.condition(end+1) = {'hit&~early&~stim.enable&autowater'};                          %  hit, no stim, aw on

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

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

meta = [];

% --- ALM --- 
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);

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
%% Null and Potent Space

clearvars -except obj meta params me sav

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -- Calculate coding directions from null and potent spaces --
% -----------------------------------------------------------------------

for sessix = 1:numel(meta)
    % -- input data (the FR of each individual neuron will be converted to
    % a value between 0 and 1)
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat, obj(sessix));

    % -- Calculate the null and potent spaces for each session
    cond2use = [2 3 7 8];   % All 2AFC hit/miss trials, all AW hit/miss trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
    nullalltime = 0;        % use all time points to estimate null space if 1
    AWonly = 0;             % use only AW to find null and potent spaces 
    delayOnly = 0;          % use only delay period to find null and potent spaces
    cond2proj = 2:13;       % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj,nullalltime,AWonly,delayOnly);
end
%%
modparams.quals2excl = {'Poor','Multi','Noisy'};
modparams.sm = 30;                                                    % Amount that you want to smooth PSTHs by
conds2plot = [12,13];
psth_all = [];
for sessix = 1:length(meta)
    sessCells = obj(sessix).psth(:,:,conds2plot);
    psth_all = cat(2,psth_all,sessCells);
end
%%
% Plot PSTHs for 2AFC and AW for each cell
% for cellix = 1:size(psth_all,2)
%     for cond = 1:2
%         plot(obj(1).time,psth_all(:,cellix,cond)); hold on;
%     end
%     legend('2AFC','AW')
%     hold off;
%     pause
% end
%% Sort the cells in order of "amount of ramping"
% Take the avg firing rate for 0.1 seconds after delay onset
% Take avg firing rate for 0.1 seconds before the go cue
% Find the difference between early delay FR and late delay FR
% Cells with largest difference are first 
e1start = find(obj(1).time>-0.9,1,'first');
e1stop = e1start+40;
e2stop = find(obj(1).time<0,1,'last');
e2start = e2stop-40;

cond2use = 2;
delaydif = NaN(1,size(psth_all,2));
for cellix = 1:size(psth_all,2)
    e1avg = mean(psth_all(e1start:e1stop,cellix,cond2use),1,'omitnan');
    e2avg = mean(psth_all(e2start:e2stop,cellix,cond2use),1,'omitnan');
    delaydif(cellix) = e2avg-e1avg;
end

[~,sortix] = sort(delaydif,'descend');
sorted_psth = squeeze(psth_all(:,sortix,cond2use));
sorted_psth = mySmooth(sorted_psth,31);
%% Plot heatmap of all cells on AW trials (condition-averaged PSTH)
imagesc(obj(1).time,1:size(psth_all,2),sorted_psth')
colorbar()
colormap('jet')
xlim([-3 0])
clim([0 40])
ylabel('Neuron')
xlabel('Time from go cue (s)')
%% Find ramping cells (spk counts in early delay are significantly less than spk counts in late delay)
rampingCells = [];
cond2use = 13;
modparams.quals2excl = {'Poor','Multi','Noisy'};
modparams.sm = 30;                                                    % Amount that you want to smooth PSTHs by
modparams.measure = 'FR';                                             % Whether you want to compare firing rate ('FR') or spike counts ('spkCnt') 
modparams.subTrials = 35;
modparams.times.e1 = [e1start e1stop];  
modparams.times.e2 = [e2start e2stop]; 
modparams.sig = 0.05;
for sessix = 1:length(meta)                                 % For every session...   
    rampingcells = identifyRampingCells(sessix,obj,modparams,meta,cond2use,params);
    obj(sessix).rampingCells = rampingcells;
end
%%
function rampingcells = identifyRampingCells(sessix,obj,modparams,meta,cond2use,params)
currobj = obj(sessix);
nTrials = size(currobj.trialdat,3);
nCells = size(currobj.psth,2);
includedCells = [];
probenum = meta(sessix).probe;
spkstuff = currobj.clu{meta(sessix).probe};
epochSpikes = zeros(nCells,nTrials,2);
for c = 1:nCells                                        % For each cell...
    cellQual = currobj.clu{probenum}(c).quality;
    % Exclude cell from analysis if it is not of proper quality
    if strcmp(cellQual,modparams.quals2excl{1}) || strcmp(cellQual,modparams.quals2excl{2}) || strcmp(cellQual,modparams.quals2excl{3})
        includedCells = [includedCells,0];
    else
        for t = 1:nTrials                                                   % Go through all of the trials
            for tim = 1:length(fieldnames(modparams.times))
                spikeix = find(currobj.clu{probenum}(c).trial==t);              % Find the spikes for this cell that belong to the current trial
                spktms = currobj.clu{probenum}(c).trialtm_aligned(spikeix);     % Get the aligned times within the trial that the spikes occur
                switch tim
                    case 1
                        epoch = 'e1';
                    case 2
                        epoch = 'e2';
                end
                epochtimes = modparams.times.(epoch);
                tempspks = length(find(spktms>epochtimes(1)&spktms<epochtimes(2)));                            % Take the spikes which occur before the sample tone
                if ~isempty(tempspks)
                    epochSpikes(c,t,tim) = tempspks;                               % Save this number
                end
            end
        end
        includedCells = [includedCells,1];
    end
end

% Whether you are comparing the firing rates of cells across contexts or spike counts
if strcmp(modparams.measure,'FR')
    psth = obj(sessix).trialdat;                                 % Get the trial PSTH (time x cells x trials)
    orgpsth = NaN(nCells,nTrials,length(fieldnames(modparams.times)));  % Cells x trials x epochs   
    for tim = 1:length(fieldnames(modparams.times))
        switch tim
            case 1
                epoch = 'e1';
            case 2
                epoch = 'e2';
        end
    epochtimes = modparams.times.(epoch);
    temp = mean(psth(epochtimes(1):epochtimes(2),:,:),1);                 % Take the average FR for all cells during the presamp period
    orgpsth(:,:,tim) = temp;                                              % (cells x trials x epoch)
    end
    clear temp
    temp = orgpsth;
elseif strcmp(modparams.measure,'spkCnt')
    temp = epochSpikes;                                        % (cells x trials x epoch)
end
trix = params(sessix).trialid{cond2use};
trix2use = randsample(trix,modparams.subTrials);
for tim = 1:length(fieldnames(modparams.times))
    epochAvg{tim} = temp(:,trix2use,tim);
end

% The p-value that you want to perform the ranksum test at
[hyp] = getContextModulatedCells(epochAvg,modparams.sig);

rampingcells = find(hyp&includedCells);
clear epochAvg
clear includedCells
end
