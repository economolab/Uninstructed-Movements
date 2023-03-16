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
    % -- input data (the FR of each individual neuron will be converted to
    % a value between 0 and 1)
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
%% Find average context selectivity for all dimensions in the null and potent space 
clearvars -except rez obj meta params trialstart samp

start = find(obj(1).time>trialstart,1,'first');
stop = find(obj(1).time<samp,1,'last');
conds2use = [5,10];       % With reference to the conditions which were projected onto the N/P dims above
sm = 100;                                        % Smoothing parameter for the condition-averaged projections

spacename = 'N_null_psth';
space = 'null';
[NullSel,NullNonSel,ctxtSelect_bySess.null] = findContextSelectivity_allDims(rez,meta,space,spacename,conds2use,sm,start,stop);
spacename = 'N_potent_psth';
space = 'potent';
[PotentSel,PotentNonSel,ctxtSelect_bySess.potent] = findContextSelectivity_allDims(rez,meta,space,spacename,conds2use,sm,start,stop);
%%
clearvars -except colors meta PotentSel NullSel obj params rez samp trialstart 
epochs.trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
epochs.sample = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
epochs.delay = median(obj(1).bp.ev.delay)-median(obj(1).bp.ev.(params(1).alignEvent));
epochs.go = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent));


Nullslope = findSelDimSlopes(NullSel.flipped,epochs,obj);
Potentslope = findSelDimSlopes(PotentSel.flipped,epochs,obj);
%%
colors = getColors_Updated();
X = categorical({'Presample','Delay','Response'});
X = reordercats(X,{'Presample','Delay','Response'});
Y = [mean(Nullslope.presample), mean(Potentslope.presample); ...
    mean(Nullslope.delay), mean(Potentslope.delay);...
    mean(Nullslope.resp), mean(Potentslope.resp)];
b = bar(X,Y);
b(1).FaceColor = colors.null;
b(2).FaceColor = colors.potent;
hold on; 
ylabel('Slope (a.u)')
xlabel('Epoch')
legend('Null','Potent','Location','best')

%% Functions
function [allSessSel,allSessNonSel,ctxtSelect_bySess] = findContextSelectivity_allDims(rez,meta,space,spacename,conds2use,sm,start,stop)

ctxtSelect_bySess = NaN(size(rez(1).N_potent,1),length(meta));
allSessSel.flipped = [];
allSessSel.nonflip = [];
allSessNonSel = [];
for sessix = 1:numel(meta)                      % For each session...
    % For null and potent spaces

    % Get the trial-averaged projections for all conditions onto each null or potent dimension
    proj = rez(sessix).(spacename);
    nDims = size(proj,2);

    % Take the difference in projections for 2AFC and AW (2AFC-AW)
    selAllDims = mySmooth(proj(:,:,conds2use(1)),sm)-mySmooth(proj(:,:,conds2use(2)),sm);
    tempSel = [];

    % Only consider selectivity for selective dimensions; Flip
    % dimensions which have negative selectivity in the presample period
    for d = 1:nDims                                                 % For every dimension...
        if rez(sessix).selectiveDimensions.(space)(d)               % If it is a selective dimension...
            temp = selAllDims(:,d);                                 % Take the selectivty for that dimension
            preSamptemp = mean(temp(start:stop,:),1);               % Find the average presample selectivity for this dimension
            allSessSel.nonflip = [allSessSel.nonflip,temp];
            if preSamptemp<0                                        % If the selectivity is negative in the presample period
                temp = -1*(temp);                                   % Flip it
            end
            tempSel = [tempSel,temp];
            allSessSel.flipped = [allSessSel.flipped,temp];                         % Store the selectivity in all dimensions, not averaged by session
        else
            temp = selAllDims(:,d);                                 % Store the selectivity in non-selective dimensions
            allSessNonSel = [allSessNonSel,temp];
        end
    end
    ctxtSelect_bySess(:,sessix) = mean(tempSel,2,'omitnan');       % Store the average selectivity across all dimensions in a session
end
end



%% Plotting functions
