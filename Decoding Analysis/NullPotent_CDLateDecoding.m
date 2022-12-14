% DECODING Null and potent CDlate FROM ALL KINEMATIC FEATURES
clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Decoding Analysis'));
%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 0; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off
params.condition(end+1) = {'R&no&~stim.enable&~autowater'};              % no right, no stim, aw off
params.condition(end+1) = {'L&no&~stim.enable&~autowater'};              % no left, no stim, aw off

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = (1/100)*3;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
% params.quality = {'Excellent','Great','Good','Fair','Multi'};

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

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
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat);

    % -- null and potent spaces
    message = strcat('----Calculating N/P spaces for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    cond2use = [2 3 4 5]; % right hit, left hit, right miss, left miss
    cond2proj = [2:7,1];  % All conditions, last entry = all single trials
    nullalltime = 0; % use all time points to estimate null space if 1
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj, nullalltime);
    
    % -- coding dimensions
    message = strcat('----Calculating N/P CDs for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    cond2use = [1 2]; % right hits, left hits (corresponding to null/potent psths in rez)
    cond2proj = [1:6]; % right hits, left hits, right miss, left miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
    cd_null(sessix) = getCodingDimensions(rez(sessix).N_null_psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);
    cd_potent(sessix) = getCodingDimensions(rez(sessix).N_potent_psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);

end

cd_null_all = concatRezAcrossSessions(cd_null);
cd_potent_all = concatRezAcrossSessions(cd_potent);
%% Load kinematic data for each session

for sessix = 1:nSessions
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% Project single trials onto Null and Potent CDs
disp('----Projecting single trials onto CDlate----')
cd = 'late';

[cd_null,cd_potent] = getNPSingleTrialProjs(obj,cd,cd_null,cd_potent,rez);
%% Predict Null and Potent CDLate from all DLC and ME features
% params
regr.nFolds = 4;                                     % number of iterations (bootstrap)

regr.binSize = 30;                                   % Bin size to decode over (in milliseconds)
difTime = params(1).dt*1000;                        % Convert time-difference from s to ms
regr.dt = floor(regr.binSize / difTime);              % How many samples confer your bin size
regr.tm = obj(1).time(1:regr.dt:numel(obj(1).time));  % Create a new time axis over which to decode (based on the bin size that you want)
regr.numT = numel(regr.tm);                           % Number of time-points

regr.train = 1;                                      % fraction of trials to use for training (1-train for testing)

% match number of right and left hits, and right and left misses
cond2use = 2:5;
hitcond = [1 2];
misscond = [3 4];

mask = true(numel(kin(1).featLeg),1);
regr.featix = find(mask);                    % Indices in the feature legend that correspond to specified body part

% Do the decoding for each session
for sessix = 1:numel(obj)
    disp(['Decoding session ' num2str(sessix) ' / ' num2str(numel(obj))])
    %%% NULL %%%
    string = 'null';
    space = cd_null;
    [decodedNull(sessix).trueVals,decodedNull(sessix).modelpred] = doCDLateDecoding(params,sessix,cond2use,hitcond,misscond,space,kin,regr,string);

    %%% POTENT %%%
    string = 'potent';
    space = cd_potent;
    [decodedPotent(sessix).trueVals,decodedPotent(sessix).modelpred] = doCDLateDecoding(params,sessix,cond2use,hitcond,misscond,space,kin,regr,string);
end
%% Find and plot the correlation of predicted and true Null/Potent CDlate
delR2.null = [];
delR2.potent = [];
% Each dot = an average value of CDlate during the delay period 
start = find(obj(1).time>-0.9,1,'first');
stop = find(obj(1).time<-0.05,1,'last');

for sessix = 1:length(meta)
    string = 'null';
    decoded = decodedNull;
    figure()
    R2 = plotNP_PredictedCDLate_Scatter(decoded, sessix, meta, string,start,stop);
    delR2.null = [delR2.null, R2];

    string = 'potent';
    decoded = decodedPotent;
    R2 = plotNP_PredictedCDLate_Scatter(decoded, sessix, meta, string,start,stop);
    delR2.potent = [delR2.potent, R2];
end
%% Scatter plot of R2 values on each session between true and predicted Null and Potent CDlate
figure();
scatter(delR2.null, delR2.potent)
xlabel('R2, true and pred Null CDlate')
ylabel('R2, true and pred Potent CDlate')
%% Bar plot of R2 values on each session between true and predicted Null and Potent CDlate
figure();
X = categorical({'Null','Potent'});
X = reordercats(X,{'Null','Potent'});
Y = [mean(delR2.null), mean(delR2.potent)];
b = bar(X,Y); hold on;
b.FaceColor = 'flat';
b.CData(1,:) = [.1 0.6 .1];
b.CData(2,:) = [1 0.2 1];
scatter(1,delR2.null,35,[0.55 0.55 0.55],'filled','MarkerEdgeColor','black')
scatter(2,delR2.potent,35,[0.55 0.55 0.55],'filled','MarkerEdgeColor','black')
ylabel("R^2 value")
%% Plot an example session of CDlate prediction vs true value
for sessix = 11%1:length(meta)
    sesstitle = strcat(meta(sessix).anm, {' '},meta(sessix).date);

    

    % Plot average CDlate projections and predictions for the session
    alpha = 0.2;
    
    % Calculate averages and standard deviation for true CD and predicted CD
    % for this session
    [avgCD,stdCD] = getAvgStd(decodedNull(sessix));
    subplot(1,2,1)
    plotExampleCDLatePrediction(obj,regr,avgCD,stdCD,sesstitle,delR2.null(sessix),alpha,decodedNull(sessix))

    % Calculate averages and standard deviation for true CD and predicted CD
    % for this session
    [avgCD,stdCD] = getAvgStd(decodedPotent(sessix));
    subplot(1,2,2)
    plotExampleCDLatePrediction(obj,regr,avgCD,stdCD,sesstitle,delR2.potent(sessix),alpha,decodedPotent(sessix))
end
%% FUNCTIONS
function [avgCD,stdCD] = getAvgStd(decodedNull)
avgCD.Rhit.true = mean(decodedNull.trueVals.Rhit,2,'omitnan');      % Get average true CDlate for R and L hits for this session
avgCD.Lhit.true = mean(decodedNull.trueVals.Lhit,2,'omitnan');
stdCD.Rhit.true = std(decodedNull.trueVals.Rhit,0,2,'omitnan');     % Get standard deviation of true CDlate for R and L hits
stdCD.Lhit.true = std(decodedNull.trueVals.Lhit,0,2,'omitnan');

decodedNull.modelpred.Rhit = fillmissing(decodedNull.modelpred.Rhit,"nearest");
decodedNull.modelpred.Lhit = fillmissing(decodedNull.modelpred.Lhit,"nearest");
infix = find(isinf(decodedNull.modelpred.Rhit)); decodedNull.modelpred.Rhit(infix) = 0;
infix = find(isinf(decodedNull.modelpred.Lhit)); decodedNull.modelpred.Lhit(infix) = 0;
avgCD.Rhit.pred = mean(decodedNull.modelpred.Rhit,2,'omitnan');     % Get average predicted CDlate for R and L hits for this session
avgCD.Lhit.pred = mean(decodedNull.modelpred.Lhit,2,'omitnan');
stdCD.Rhit.pred = std(decodedNull.modelpred.Rhit,0,2,'omitnan');    % Get stdev of predicted CDlate for R and L hits for this session
stdCD.Lhit.pred = std(decodedNull.modelpred.Lhit,0,2,'omitnan');
end

function plotExampleCDLatePrediction(obj,regr,avgCD,stdCD,sesstitle,R2,alpha,decodedNull)
plot(obj(1).time,avgCD.Rhit.true,'Color','blue','LineWidth',2); hold on; ax = gca;
err = 1.96*(stdCD.Rhit.true/size(decodedNull.trueVals.Rhit,2));
%shadedErrorBar(obj(1).time,avgCD.Rhit.true,err,{'Color','blue','LineWidth',2}, alpha, ax)

plot(regr.tm(1:end-1),avgCD.Rhit.pred,'Color',[0.5 0.5 1],'LineStyle','--','LineWidth',2); ax = gca;
err = 1.96*(stdCD.Rhit.pred/size(decodedNull.modelpred.Rhit,2));
%shadedErrorBar(regr.tm(1:end-1),avgCD.Rhit.pred,err,{'Color',[0.5 0.5 1],'LineWidth',2}, alpha, ax)

plot(obj(1).time,avgCD.Lhit.true,'Color','red','LineWidth',2); hold on;
err = 1.96*(stdCD.Lhit.true/size(decodedNull.trueVals.Lhit,2));
%shadedErrorBar(obj(1).time,avgCD.Lhit.true,err,{'Color','red','LineWidth',2}, alpha, ax)


plot(regr.tm(1:end-1),avgCD.Lhit.pred,'Color',[1 0.5 0.5],'LineStyle','--','LineWidth',2);
err = 1.96*(stdCD.Lhit.pred/size(decodedNull.modelpred.Lhit,2));
%shadedErrorBar(regr.tm(1:end-1),avgCD.Lhit.pred,err,{'Color',[1 0.5 0.5],'LineWidth',2}, alpha, ax)


xline(0,'black','LineStyle','--','LineWidth',1.1)
xline(-0.9,'black','LineStyle','-.','LineWidth',1.1)
xline(-2.5,'black','LineStyle','-.','LineWidth',1.1)

legend('R hit true','R hit predicted','L hit true','L hit predicted','Location','best')
ylabel('a.u.')
xlabel('Time since go-cue (s)')
sesstit = [sesstitle, 'R^2 = ', num2str(R2)];
title(sesstit)
xlim([-2.6 2.5])
end