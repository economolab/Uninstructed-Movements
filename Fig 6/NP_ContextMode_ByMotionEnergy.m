% Finding CDContext from neural activity that resides within the Null and Potent spaces
% Then finding selectivity between CDContext from full neural pop and
% CDContext found from null/potent reconstructions
clear,clc,close all

% add paths
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig3')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 3';
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context switching')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 6';
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context_funcs')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 5';
addpath(genpath(fullfile(figpth,'funcs')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 2';
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
%%
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
% -- Calculate coding directions from full neural population --
% -----------------------------------------------------------------------

for sessix = 1:numel(meta)
    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat, obj(sessix));

    % -- Calculate the null and potent spaces for each session
    cond2use = [2 3 4 5];   % All 2AFC hit/miss trials, all AW hit/miss trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
    nullalltime = 0;        % use all time points to estimate null space if 1
    AWonly = 0;             % use only AW to find null and potent spaces 
    delayOnly = 0;          % use only delay period to find null and potent spaces
    cond2proj = [6 7];      % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj,nullalltime,AWonly,delayOnly);

    % -- Find coding dimensions from RECONSTRUCTED full neural activity which is reconstructed from the null and potent spaces
    cond2use = [1 2];            % (NUMBERING ACCORDING TO THE CONDITIONS PROJECTED INTO NULL AND POTENT SPACES, i.e. which of the conditions specified in 'cond2proj' above do you want to use?)
    cond2proj = [1 2];           % 2AFC hits/misses, AW hits/misses(corresponding to null/potent psths in rez)
    cond2use_trialdat = [6 7];   % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    cd_null(sessix) = getCodingDimensions_Context(rez(sessix).recon_psth.null,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);
    cd_potent(sessix) = getCodingDimensions_Context(rez(sessix).recon_psth.potent,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);

    % Calc CDContext from full neural pop
    cd_context(sessix) = getCodingDimensions_Context(obj(sessix).psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);
end
%% Get single trial projections onto CDContext
orthogonalize = 'non-orthog';                                       % Set to orthogonalize if you want the projections to be onto the orthogonalized CDs
disp('----Projecting single trials onto CDContext----')
cd = 'context';
cd_context = getSingleTrialProjs(cd_context,obj,cd,orthogonalize);
%% Get single trial NP projections onto CDContext
cd = 'context';
[cd_null,cd_potent] = getNPSingleTrialProjs(obj,cd,cd_null,cd_potent,rez);
%% Get trials with a lot of presample motion energy and little presample motion energy
clearvars -except obj meta params me sav rez cd_null cd_potent cd_context

% Find the times corresponding to trial start and the sample period
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
startix = find(obj(1).time>trialstart,1,'first');
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
stopix = find(obj(1).time<samp,1,'last');

trials2cutoff = 40;                 % Trials to discount at the end of the session (for motion energy)
cond2use = [6,7];                   % 2AFC trials and AW trials
ngroups = 3;
for sessix = 1:length(meta)
    nTrials = obj(sessix).bp.Ntrials;
    cutoff = nTrials-trials2cutoff;
    for c = 1:length(cond2use)
        if c==1
            trialContext = 'afc';
        else
            trialContext = 'aw';
        end
        cond = cond2use(c);
        % Get trials for current cond
        trix2use = params(sessix).trialid{cond};
        trix2use = trix2use(trix2use<cutoff);           % Only use trials that come before the 'end of session cutoff'

        % Get ME and projections onto Ramping Mode for these trials
        MEtrix = me(sessix).data(:,trix2use);
        context = cd_context(sessix).singleProj(:,trix2use);

        % Find avg ME during presamp on each trial
        avgME = mean(MEtrix(startix:stopix,:),1,'omitnan');
        % Sort ME and Context proj by ME descending order
        [~,sortix] = sort(avgME,'descend');
        sortedME = MEtrix(:,sortix);
        sortedcontext = context(:,sortix);

        nTrials = length(sortix);
        trixPerGroup = floor(nTrials/ngroups);                  % How many trials you want to be in each group
        cnt = 1;
        tempME = cell(1,ngroups);
        tempCont = cell(1,ngroups);
        for g = 1:ngroups
            if g==ngroups
                ixrange = cnt:nTrials;
            else
                ixrange = cnt:(cnt+trixPerGroup);
            end
            tempME{g} = mean(sortedME(:,ixrange),2,'omitnan');
            tempCont{g} = mean(sortedcontext(:,ixrange),2,'omitnan');
            cnt = cnt+trixPerGroup+1;
        end
    end
    grouped(sessix).ME.(trialContext) = tempME;
    grouped(sessix).cont.(trialContext) = tempCont;
end
%%

%% Plotting functions
function plotAvgSelectivity_NP(NullSel,PotentSel, meta, colors, alpha, smooth, trialstart, samp,obj,samples)
figure();
ax = gca;
col = colors.null;
% temperr = 1.96*(mySmooth(std(ctxtSelect.null,0,2),smooth)/sqrt(length(meta)));
% toplot = mySmooth(mean(ctxtSelect.null,2,'omitnan'),smooth);
% What to normalize by for the confidence interval calculations
if strcmp(samples,'sessions')
    normaliz = length(meta);
elseif strcmp(samples,'alldims')
    normaliz = size(NullSel.flipped,2);
end
%temperr = 1.96*(std(NullSel.flipped,0,2)/sqrt(normaliz));
temperr = (std(NullSel.flipped,0,2)/sqrt(normaliz));    % SEM
toplot = mean(NullSel.flipped,2,'omitnan');
shadedErrorBar(obj(1).time,toplot,temperr,{'Color',col,'LineWidth',2}, alpha, ax); hold on;


if strcmp(samples,'sessions')
    normaliz = length(meta);
elseif strcmp(samples,'alldims')
    normaliz = size(PotentSel.flipped,2);
end

col = colors.potent;
%temperr = 1.96*(std(PotentSel.flipped,0,2)/sqrt(normaliz));
temperr = (std(PotentSel.flipped,0,2)/sqrt(normaliz));      % SEM
toplot = mean(PotentSel.flipped,2,'omitnan');
shadedErrorBar(obj(1).time,toplot,temperr,{'Color',col,'LineWidth',2}, alpha, ax); hold on;

xlim([trialstart 2])
xline(samp,'k--','LineWidth',1)
xline(samp+1.3,'k--','LineWidth',1)
xline(samp+1.3+0.9,'k--','LineWidth',1)
xlabel('Time from go cue/water drop (s)')
ylabel('Selectivity across all null/potent dimensions (a.u.)')
end

function plotSelectivityHeatmap(delay,go,samp,selectNorm,obj,ContextColorMap)
subplot(1,2,1)
nDims = size(selectNorm.null,2);
imagesc(obj(1).time,1:nDims,selectNorm.null'); c = colorbar;
ax = gca;
colormap(ax,linspecer)
title('Null')
ylabel('Dimensions')
xlabel('Time from go cue (s)')
line([samp,samp],[1,nDims],'Color','black','LineStyle','--','LineWidth',1.5)
line([delay,delay],[1,nDims],'Color','black','LineStyle','--','LineWidth',1.5)
line([go,go],[1,nDims],'Color','black','LineStyle','--','LineWidth',1.5)

subplot(1,2,2)
nDims = size(selectNorm.potent,2);
imagesc(obj(1).time,1:nDims,selectNorm.potent'); c = colorbar;
ylabel(c,'Normalized selectivity','FontSize',9,'Rotation',90);
ax = gca;
colormap(ax,linspecer)
title('Potent')
line([samp,samp],[1,nDims],'Color','black','LineStyle','--','LineWidth',1.5)
line([delay,delay],[1,nDims],'Color','black','LineStyle','--','LineWidth',1.5)
line([go,go],[1,nDims],'Color','black','LineStyle','--','LineWidth',1.5)
xlabel('Time from go cue (s)')
end

