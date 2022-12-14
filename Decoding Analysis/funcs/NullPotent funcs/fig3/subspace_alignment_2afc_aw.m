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

clc

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

params.condition(end+1) = {'R&hit&~stim.enable&autowater'};             % right hits, no stim, aw on
params.condition(end+1) = {'L&hit&~stim.enable&autowater'};             % left hits, no stim, aw on
params.condition(end+1) = {'R&miss&~stim.enable&autowater'};            % error right, no stim, aw on
params.condition(end+1) = {'L&miss&~stim.enable&autowater'};            % error left, no stim, aw on


params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

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

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];


% autowater sessions—it’s all the sessions from JEB6, JEB7, EKH1, EKH3, JGR2, JGR3
% --- ALM --- 
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);
% meta = loadJEB15_ALMVideo(meta,datapth);


% --- M1TJ ---
% meta = loadJEB13_M1TJVideo(meta,datapth);

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

clearvars -except obj meta params me sav datapth

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -- subspace alignment b/w contexts --
% -----------------------------------------------------------------------

for sessix = 1:numel(meta)
    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat);

    % -- null and potent spaces
    cond2use = [2 3 4 5]; % right hit, left hit, right miss, left miss (afc)
    cond2proj = [2:3];
    nullalltime = 0; % use all time points to estimate null space if 1
    rez_afc(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj, nullalltime);

    cond2use = [6 7 8 9]; % right hit, left hit, right miss, left miss (aw)
    cond2proj = [6:7];
    nullalltime = 0; % use all time points to estimate null space if 1
    rez_aw(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj, nullalltime);

    % subspace alignment between contexts
    % This procedure reveals the amount of quiet-epoch and movement-epoch
    % variance is shared between contexts
    Q.null.aw = rez_aw(sessix).Qnull;
    Q.null.afc = rez_afc(sessix).Qpotent;
    Q.potent.aw = rez_aw(sessix).Qpotent;
    Q.potent.afc = rez_afc(sessix).Qpotent;
    C.null.aw = rez_aw(sessix).covNull;
    C.null.afc = rez_afc(sessix).covNull;
    C.potent.aw = rez_aw(sessix).covPotent;
    C.potent.afc = rez_afc(sessix).covPotent;
    sigma.null.aw = sort(eig(C.null.aw),'descend'); sigma.null.aw = sigma.null.aw(1:rez_aw(sessix).dPrep);
    sigma.null.afc = sort(eig(C.null.afc),'descend'); sigma.null.afc = sigma.null.afc(1:rez_afc(sessix).dPrep);
    sigma.potent.aw = sort(eig(C.potent.aw),'descend'); sigma.potent.aw = sigma.potent.aw(1:rez_aw(sessix).dMove);
    sigma.potent.afc = sort(eig(C.potent.afc),'descend'); sigma.potent.afc = sigma.potent.afc(1:rez_afc(sessix).dMove);
    ai.null.afc_aw(sessix) = get_alignment_index(Q.null.afc,C.null.aw,sigma.null.aw);
    ai.null.aw_afc(sessix) = get_alignment_index(Q.null.aw,C.null.afc,sigma.null.afc);
    ai.potent.afc_aw(sessix) = get_alignment_index(Q.potent.afc,C.potent.aw,sigma.potent.aw);
    ai.potent.aw_afc(sessix) = get_alignment_index(Q.potent.aw,C.potent.afc,sigma.potent.afc);

    % subspace alignment between epochs
    epochai(sessix) = get_alignment_index(Q.potent.afc,C.null.afc,sigma.null.afc);
    
% %     ang.null(sessix) = rad2deg(subspace(Q.null.afc,Q.null.aw));
% %     ang.potent(sessix) = rad2deg(subspace(Q.potent.afc,Q.potent.aw));

    % cosine similarity b/w cov matrices (1 equal, 0 orthogonal, -1
    % antiequal)
%     a = C.null.afc(:);
%     b = C.potent.afc(:);
%     cosSim(sessix) = dot(a,b)/(norm(a)*norm(b));


%     % correlation matrix distance 
%     corrnull = corr(rez_afc(sessix).N.null);
%     corrpotent = corr(rez_afc(sessix).N.potent);
%     cmd2(sessix) = 1 - trace(corrnull * corrpotent) ./ ( norm(corrnull,'fro') * norm(corrpotent,'fro')  );
%     cmd(sessix) = 1 - trace(rez_afc(sessix).covNull * rez_afc(sessix).covPotent) ./ ...
%           ( norm(rez_afc(sessix).covNull,'fro') * norm(rez_afc(sessix).covPotent,'fro') );
end


%% plots

close all

% % - ve
plotVarianceExplained_NP(rez_afc);

plotSubspaceAlignment(ai);
% plotCosSim(cosSim);



cols = [0.1 0.1 0.1; 0.4 0.4 0.4];
f=figure; hold on;
ax = gca;
div = 1;
xs = [1];
i = 1; % potent afc_aw 
temp = epochai;
h(i) = bar(xs(i), mean(temp));
h(i).FaceColor = cols(i,:);
h(i).EdgeColor = 'none';
h(i).FaceAlpha = 0.5;
scatter(xs(i)*ones(size(temp)),temp,60,'MarkerFaceColor',cols(i,:)./div, ...
        'MarkerEdgeColor','k','LineWidth',1,'XJitter','randn','XJitterWidth',0.25, ...
        'MarkerFaceAlpha',0.7)
errorbar(h(i).XEndPoints,mean(temp),std(temp)./sqrt(numel(temp)),'LineStyle','none','Color','k','LineWidth',1);
ylim([0 1])


















