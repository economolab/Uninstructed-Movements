% Finding CDContext from neural activity that resides within the Null and Potent spaces
% Then finding selectivity between CDContext from full neural pop and
% CDContext found from null/potent reconstructions
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
addpath(genpath(fullfile(utilspth,'figNP')));
figpth = [basepth  '\Uninstructed-Movements\Fig 3'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context switching')));
figpth = [basepth  '\Uninstructed-Movements\Fig 6'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context_funcs')));
figpth = [basepth  '\Uninstructed-Movements\Fig 5'];
addpath(genpath(fullfile(figpth,'funcs')));
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
params.condition(end+1) = {'hit&~stim.enable&~autowater'};               % all 2AFC hits, no stim
params.condition(end+1) = {'hit&~stim.enable&autowater'};                % all AW hits, no stim
params.condition(end+1) = {'miss&~stim.enable&~autowater'};              % error 2AFC, no stim, aw off
params.condition(end+1) = {'miss&~stim.enable&autowater'};               % error AW, no stim

params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % all 2AFC hits, ~early, no stim
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};                % all AW hits, ~early,no stim

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

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
meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
%meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);
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
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end
%% Find Null and Potent Spaces using all trials (no train/test split yet)
clearvars -except obj meta params me sav

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -----------------------------------------------------------------------
for sessix = 1:numel(meta)
    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix));
    zscored(sessix).trialdat =  trialdat_zscored;

    % -- Calculate the null and potent spaces for each session
    cond2use = [2 3 4 5];   % All 2AFC hit/miss trials, all AW hit/miss trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
    nullalltime = 0;        % use all time points to estimate null space if 1
    AWonly = 0;             % use only AW to find null and potent spaces 
    delayOnly = 0;          % use only delay period to find null and potent spaces
    cond2proj = [6 7];      % 2AFC hits (non-early); AW hits (non-early)
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj,nullalltime,AWonly,delayOnly); 
end
%% Sum of the squared magnitude across dimensions 
% rez.N_space = (time x trials x dims)
% Want to get the sum-squared magnitude within the null and potent spaces
% sumsq_mag.space = (time x trials) where each entry is the sum of the squared magnitude across dimensions
for sessix = 1:length(meta)
    sq_mag_null = rez(sessix).N_null.^2;                                    % Square each element (stays time x trials x dims) -- don't care about the polarity
    sumsq_mag(sessix).null = zscore(sum(sq_mag_null,3,'omitnan'));          % Add the squared magnitudes across all dimensions (time x trials) and zscore

    sq_mag_potent = rez(sessix).N_potent.^2; 
    sumsq_mag(sessix).potent = zscore(sum(sq_mag_potent,3,'omitnan'));
end
%% SANITY CHECK -- cross-correlation between motion energy and null or potent space 
sm = 21;

dt = params(1).dt;
maxlag = int64(3 ./ dt);                                                                       % Specify which lags you want to compute the cross-correlation for

laglength = (2*maxlag)+1;
allcorrsME.null = NaN(laglength,length(meta));
allcorrsME.potent = NaN(laglength,length(meta));
for sessix = 1:length(meta)
    nTrials = size(sumsq_mag(sessix).null,2);
    for trix = 1:nTrials
        currME = zscore(me(sessix).data(:,trix));                                              % z-score the motion energy (ask Munib why this matters?)
        currPot = sumsq_mag(sessix).potent(:,trix);
        currNull = sumsq_mag(sessix).null(:,trix);

%         currME = mySmooth(zscore(me(sessix).data(:,trix)),sm);                                              % z-score the motion energy (ask Munib why this matters?)
%         currPot = mySmooth(sumsq_mag(sessix).potent(:,trix),sm);
%         currNull = mySmooth(sumsq_mag(sessix).null(:,trix),sm);

        [Ncorrs(:,trix),lagtm] = xcorr(currME,currNull,maxlag,'normalized'); % {session}(time,trial,dim)
        [Pcorrs(:,trix),lagtm] = xcorr(currME,currPot,maxlag,'normalized'); % {session}(time,trial,dim)
        % 'Normalized' parameter normalizes the correlation values so that the auto-correlations at time-zero = 1

        %%% Sanity check -- what do the cross-correlations between ME and
        %%% null or ME and potent look like on single trials
        %         subplot(2,1,1)
        %         plot(currNull); hold on; plot(currPot);plot(currME,'LineWidth',1.5); legend('Null','Potent','ME'); hold off
        %         subplot(2,1,2)
        %         plot(lagtm, Ncorrs(:,trix)); hold on; plot(lagtm,Pcorrs(:,trix)); legend('Null','Potent'); hold off
        %         pause
    end
    allcorrsME.null(:,sessix) = mean(Ncorrs,2,'omitnan');
    allcorrsME.potent(:,sessix) = mean(Pcorrs,2,'omitnan');
end
%% Find the cross-correlation between null and potent on single trials
dt = params(1).dt;
maxlag = int64(3 ./ dt);                                                                       % Specify which lags you want to compute the cross-correlation for

laglength = (2*maxlag)+1;
allcorrsNP.true = NaN(laglength,length(meta));
allcorrsNP.shuff = NaN(laglength,length(meta));
for sessix = 1:length(meta)
    nTrials = size(sumsq_mag(sessix).null,2);
    randtrials = randsample(1:nTrials,nTrials,'false');                                       % Shuffle the trial order as a control
    for trix = 1:nTrials
        shufftrix = randtrials(trix);                                                         % Get the shuffled trial number
        
        currPot = sumsq_mag(sessix).potent(:,trix);                                           % Get the potent space proj for the current trial
        currNull = sumsq_mag(sessix).null(:,trix);                                            % Null space proj for the current trial
        shuffNull = sumsq_mag(sessix).null(:,shufftrix);                                      % Take the null space proj for that shuffled trial number

%         currPot = mySmooth(sumsq_mag(sessix).potent(:,trix),sm);                                           % Get the potent space proj for the current trial
%         currNull = mySmooth(sumsq_mag(sessix).null(:,trix),sm);                                            % Null space proj for the current trial
%         shuffNull = mySmooth(sumsq_mag(sessix).null(:,shufftrix),sm);                                      % Take the null space proj for that shuffled trial number


        [Tcorrs(:,trix),lagtm] = xcorr(currPot,currNull,maxlag,'normalized');                 % Cross-corr between potent and null for the same trial
        [Shuffcorrs(:,trix),lagtm] = xcorr(currPot,shuffNull,maxlag,'normalized');            % Cross-corr between potent and null on shuffled trials
        % 'Normalized' parameter normalizes the correlation values so that the auto-correlations at time-zero = 1

        %%% Sanity check -- what do the cross-correlations between null and
        %%% potent look like on single trials?
%         figure(2)
%         subplot(2,2,1)
%         plot(currNull); hold on; plot(currPot); legend('Null','Potent'); hold off
%         title('True')
%         subplot(2,2,2)
%         plot(lagtm, Tcorrs(:,trix));
%         title('True')
% 
%         subplot(2,2,3)
%         plot(shuffNull); hold on; plot(currPot); legend('ShuffNull','Potent'); hold off
%         title('Shuffled')
%         subplot(2,2,4)
%         plot(lagtm, Shuffcorrs(:,trix));
%         title('Shuffled')
%         pause
    end
    allcorrsNP.true(:,sessix) = mean(Tcorrs,2,'omitnan');                                     % Mean across trials for this session 
    allcorrsNP.shuff(:,sessix) = mean(Shuffcorrs,2,'omitnan'); 
end
%% Plot cross-corr between null or potent and ME
spacefns = {'null','potent'};
colors = getColors();
alph = 0.2;
figure(1);
subplot(1,2,1)
for ss = 1:length(spacefns)
    space = spacefns{ss};
    if strcmp(space,'null')
        col = colors.null;
    else
        col = colors.potent;
    end
    toplot = mean(allcorrsME.(space),2,'omitnan');                                              % Mean across all sessions
    err = 1.96*(std(allcorrsME.(space),0,2,'omitnan')./sqrt(length(meta)));                     % Stdev across sessions
    ax = gca;
    shadedErrorBar(lagtm*dt,toplot,err,{'Color',col,'LineWidth',2},alph,ax); hold on;
end
legend(' ',' ',' ','Null',' ',' ',' ','Potent')
xlim([-2 2])
xlabel('lag (s)')
ylabel('Correlation')
title('Cross-corr, Null/Pot and ME')
set(gca,'TickDir','out');

% Plot cross-corr btw null and potent on single trials 
shfns = {'true','shuff'};
alph = 0.2;
subplot(1,2,2)
for ss = 1:length(shfns)
    cond = shfns{ss};
    if strcmp(cond,'true')
        col = [0 0 0];
    else
        col = [0.4 0.4 0.4];
    end
    toplot = mean(allcorrsNP.(cond),2,'omitnan');
    err = 1.96*(std(allcorrsNP.(cond),0,2,'omitnan')./sqrt(length(meta)));
    ax = gca;
    shadedErrorBar(lagtm*dt,toplot,err,{'Color',col,'LineWidth',2},alph,ax); hold on;
end
legend(' ',' ',' ','Correct order',' ',' ',' ','Shuffled trials')
xlim([-2 2])
xlabel('lag (s)')
ylabel('Correlation')
title('Cross-corr, null and pot')
sgtitle('Averaged across all trials for a session; then session averaged')
ylim([-0.2 0.4])
set(gca,'TickDir','out');