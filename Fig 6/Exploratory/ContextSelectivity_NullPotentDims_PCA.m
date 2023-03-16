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
%% Concatenate condition averaged projections onto selective dimensions across sessions
conds2use = [5,10];       % With reference to the conditions which were projected onto the N/P dims above
condprojAllSess.null = [];
condprojAllSess.potent = [];
for sessix = 1:length(meta)
    space = 'null';
    spacename = 'N_null_psth';
    selDimIx = find(rez(sessix).selectiveDimensions.(space));
    temp = rez(sessix).(spacename)(:,selDimIx,conds2use);
    condprojAllSess.(space) = cat(2,condprojAllSess.(space),temp);

    space = 'potent';
    spacename = 'N_potent_psth';
    selDimIx = find(rez(sessix).selectiveDimensions.(space));
    temp = rez(sessix).(spacename)(:,selDimIx,conds2use);
    condprojAllSess.(space) = cat(2,condprojAllSess.(space),temp);
end
%% Do PCA on selective dimensions
space = 'null';
[coeff.null,score.null,latents.null,VE.null] = doPCAonSelDims(condprojAllSess,space);

space = 'potent';
[coeff.potent,score.potent,latents.potent,VE.potent] = doPCAonSelDims(condprojAllSess,space);
%% Plot projections onto top dimensions for null and potent space
dims2plot = 1:3;
colors = getColors_Updated();
figure();
plotTopPCs(dims2plot,latents,colors,VE,obj)

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

function selectNorm = normalizeSelDims(spaceSel,start, stop)
% Normalize each dimension to its own max selectivity (to show dynamics of
% selectivity for each dimension over the course of the trial)
maxSelect =  max(abs(spaceSel.nonflip),[],1);                    % Max magnitude of context selectivity for each dimension
selectNorm = [];
for c = 1:length(maxSelect)                                        % For every dimension... 
    selectNorm(:,c) = spaceSel.nonflip(:,c)./maxSelect(c);       % Normalize all selectivity values to the max selectivity val
end 
% response = mean(selectNorm.null(start:stop,:),1);
% [~,sortix] = sort(response,'descend');
full = mean(selectNorm,1);
[~,sortix] = sort(full,'descend');
selectNorm = selectNorm(:,sortix);
end

function [coeff,score,latents,VE] = doPCAonSelDims(condprojAllSess,space)
% concatenate data into a 2D array of size (C*T,D) where C is number of
% conditions, T is time points, D is number of dimensions
A = [];                                         % Pre-allocate activity matrix of dimensions
for i = 1:size(condprojAllSess.(space),3)                      % For all behavioral conditions...
    A = [A; condprojAllSess.(space)(:,:,i)];                   % Concatenate the activity matrices for all conditions 
end

Ndims = 10; % number of principal components to look at 

% coeff is the principal components. Should be of size (D,D). Each column
% is a principal component, and each element of each column is the weight
% associated with a dimension.
% score is the projection of the full dimensional data onto the principal
% components (score = psth * coeff)
[coeff, score,blah1,blah2,VE] = pca(A - mean(A), 'NumComponents', Ndims);          % Perform PCA on centered data
clear A

% project each condition onto the principal components
latents = zeros(size(condprojAllSess.(space),1),Ndims,size(condprojAllSess.(space),3));           % (time,nDims,numCond)
sm = 21;                                                                                           % How much smoothing you want to do

% The latents are the same as the 'scores' given by the PCA function
% But these are smoothed and separated in the third dimension by behavioral
% condition
for i = 1:size(condprojAllSess.(space),3)                           % For all behavioral conditions...
    for j = 1:Ndims                                                 % For all principal components...
        latents(:,j,i) = mySmooth(condprojAllSess.(space)(:,:,i)*coeff(:,j),sm);   % Assign the 'latents' to be neuron's trial-averaged PSTH weighted by their coefficients for the given PC
    end                                                             % And smooth the latent
end
end

function [allSessSelProj] = findContextProj_allDims(rez,meta,space,spacename,conds2use,sm)
allSessSelProj{1} = [];
allSessSelProj{2} = [];
for sessix = 1:numel(meta)                      % For each session...
    % For null and potent spaces

    % Get the trial-averaged projections for all conditions onto each null or potent dimension
    proj = rez(sessix).(spacename);
    nDims = size(proj,2);

    for c = 1:length(conds2use)
        condproj = mySmooth(proj(:,:,conds2use(c)),sm);

        % Only consider projections for selective dimensions
        for d = 1:nDims                                                 % For every dimension...
            if rez(sessix).selectiveDimensions.(space)(d)               % If it is a selective dimension...
                temp = condproj(:,d);                                   % Take the avg projection for that dimension for this condition

                allSessSelProj{c}= [allSessSelProj{c},temp];                  % Store the projections for all dimensions, not averaged by session
            end
        end
    end
end
end
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

function plotTopPCs(dims2plot,latents,colors,VE,obj)
for d= dims2plot
    % Null space
    subplot(2,length(dims2plot),d)
    plot(obj(1).time,latents.null(:,d,1),'Color',colors.afc,'LineWidth',1.5); hold on;
    plot(obj(1).time,latents.null(:,d,2),'Color',colors.aw,'LineWidth',1.5)
    title(['VE = ' num2str(VE.null(d)) ' %'])               % Display the variance explained by the PC as the title of the subplot
    xlim([-2.6 2.3])

    % Potent space
    subplot(2,length(dims2plot),d+3)
    plot(obj(1).time,latents.potent(:,d,1),'Color',colors.afc,'LineWidth',1.5); hold on;
    plot(obj(1).time,latents.potent(:,d,2),'Color',colors.aw,'LineWidth',1.5)
    title(['VE = ' num2str(VE.potent(d)) ' %'])             % Display the variance explained by the PC as the title of the subplot
    xlim([-2.6 2.3])
end
end