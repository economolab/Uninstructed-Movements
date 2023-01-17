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
    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat, obj(sessix));

    % -- Calculate the null and potent spaces for each session
    cond2use = [2 3 7 8];    % All 2AFC hit trials, all AW hit trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
    nullalltime = 0;      % use all time points to estimate null space if 1
    cond2proj = 2:11;     % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj,nullalltime);
end
%% Find dimensions which are selective for context
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
modparams.start = find(obj(1).time>trialstart,1,'first');
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
modparams.stop = find(obj(1).time<samp,1,'last');

modparams.subTrials = 40;
modparams.sig = 0.05;
conds2use = [6,11];       % With reference to params.trialid

for sessix = 1:length(meta)
    currez = rez(sessix);
    spacename = 'N_null';
    space = 'null';
    [selectiveDims,projdif] = findAllSelectiveDims(currez,spacename,params(sessix),modparams,conds2use);
    rez(sessix).selectiveDimensions.(space) = selectiveDims;
    rez(sessix).projdif.(space) = projdif; 

    spacename = 'N_potent';
    space = 'potent';
    [selectiveDims,projdif] = findAllSelectiveDims(currez,spacename,params(sessix),modparams,conds2use);
    rez(sessix).selectiveDimensions.(space) = selectiveDims;
    rez(sessix).projdif.(space) = projdif; 
end
%%
for sessix = 1:length(meta)
    space = 'null';
    selDims = rez(sessix).selectiveDimensions.(space);
    
end
%% Find context selectivity for all dimensions in the null and potent space 
conds2use = [5,10];       % With reference to the conditions which were projected onto the N/P dims above
ctxtSelect.null = NaN(size(rez(1).N_potent,1),length(meta));
ctxtSelect.potent = NaN(size(rez(1).N_potent,1),length(meta));

for sessix = 1:numel(meta)
    for i=1:2
        if i==1
            spacename = 'N_null_psth';
            space = 'null';
        else
            spacename = 'N_potent_psth';
            space = 'potent';
        end
        % Get the projections for all conditions onto each null or potent dimension
        proj = rez(sessix).(spacename);
        nDims = size(proj,2);

        % Take the difference projections for 2AFC and AW
        selAllDims = proj(:,:,conds2use(1))-proj(:,:,conds2use(2));
        tempSel = [];
        
        % Only consider selectivity for selective dimensions
        for d = 1:nDims                                                 % For every dimension...
            if rez(sessix).selectiveDimensions.(space)(d)                       % If it is a selective dimension...
                temp = selAllDims(:,d);                                    % Take the selectivty for that dimension
                if ~rez(sessix).projdif.(space)(d)                              % If it is selective for AW
                    temp = -1*(temp);                                   % Flip it
                end
                tempSel = [tempSel,temp];
            end
        end
        ctxtSelect.(space)(:,sessix) = mean(tempSel,2,'omitnan');  
    end
end
%%
colors = getColors_Updated();
alpha = 0.2;
figure();
ax = gca;
smooth = 31;

col = colors.null;
% temperr = 1.96*(mySmooth(std(ctxtSelect.null,0,2),smooth)/sqrt(length(meta)));
% toplot = mySmooth(mean(ctxtSelect.null,2,'omitnan'),smooth);
temperr = 1.96*(std(ctxtSelect.null,0,2)/sqrt(length(meta)));
toplot = mean(ctxtSelect.null,2,'omitnan');
shadedErrorBar(obj(1).time,toplot,temperr,{'Color',col,'LineWidth',2}, alpha, ax); hold on;

col = colors.potent;
temperr = 1.96*(std(ctxtSelect.potent,0,2)/sqrt(length(meta)));
toplot = mean(ctxtSelect.potent,2,'omitnan');
shadedErrorBar(obj(1).time,toplot,temperr,{'Color',col,'LineWidth',2}, alpha, ax); hold on;

xlim([trialstart 2.5])
xline(samp,'k--','LineWidth',1.5)
xline(samp+1.3,'k--','LineWidth',1.5)
xline(samp+1.3+0.9,'k--','LineWidth',1.5)
xlabel('Time from go cue (s)')
ylabel('Selectivity across all null/potent dimensions (a.u.)')


