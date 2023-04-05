% Find context selectivity across only Null and potent dimensions which are selective for context
% Take the difference between contexts for each selective dimension.  Take the mean
% across all selective dimensions within a session
% -------------------------------------------------------------------------------------
% Using all 2AFC and all AW trials to find the Null and Potent Spaces
% -------------------------------------------------------------------------------------
clear,clc,close all

whichcomp = 'Laptop';                                                % LabPC or Laptop

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
addpath(genpath(fullfile(utilspth,'fig3')));
figpth = [basepth '\Uninstructed-Movements\Fig 4'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context switching')));
figpth = [basepth '\Uninstructed-Movements\Fig 7'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context_funcs')));
figpth = [basepth '\Uninstructed-Movements\Fig 6'];
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

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

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
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix));

    % -- Calculate the null and potent spaces for each session
    cond2use = [2 3 7 8];   % All 2AFC hit trials, all AW hit trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
    nullalltime = 0;        % use all time points to estimate null space if 1
    AWonly = 0;             % use only AW to find null and potent spaces 
    delayOnly = 0;          % use only delay period to find null and potent spaces
    cond2proj = 2:11;       % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj,nullalltime,AWonly,delayOnly);
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
sm = 31;                                        % Smoothing parameter for the condition-averaged projections

spacename = 'N_null_psth';
space = 'null';
[NullSel,NullNonSel,ctxtSelect_bySess.null] = findContextSelectivity_allDims(rez,meta,space,spacename,conds2use,sm,start,stop);
spacename = 'N_potent_psth';
space = 'potent';
[PotentSel,PotentNonSel,ctxtSelect_bySess.potent] = findContextSelectivity_allDims(rez,meta,space,spacename,conds2use,sm,start,stop);
%% Normalize selective dimensions
start = find(obj(1).time>trialstart,1,'first');
stop = find(obj(1).time<samp,1,'last');

% Sort dimensions according to average selectivity in the presample period and normalize each dimension to its max selectivity
%%% Null space %%%
selectNorm.null = normalizeSelDims(NullSel,start, stop);
nonselectNorm.null = normalizeSelDims(NullNonSel,start, stop);

%%% Potent space %%%
selectNorm.potent = normalizeSelDims(PotentSel,start, stop);
nonselectNorm.potent = normalizeSelDims(PotentNonSel,start, stop);
%% Arrange dimensions in the order that you want to plot them
start = find(obj(1).time>trialstart,1,'first');
stop = find(obj(1).time<samp,1,'last');

Sel2Use = selectNorm.null;
NonSel2Use = nonselectNorm.null;
[toplot.Null,NullNums] = sortDimsforPlotting(Sel2Use, NonSel2Use, start,stop);

Sel2Use = selectNorm.potent;
NonSel2Use = nonselectNorm.potent;
[toplot.Potent,PotentNums] = sortDimsforPlotting(Sel2Use, NonSel2Use, start,stop);
%%
load('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\ContextColormap.mat');
tLim.start = -3;
tLim.stop = 1;
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));

nDims = size(toplot.Null,2);
disp(num2str(nDims))
nums = NullNums;
subplot(1,2,1)
imagesc(obj(1).time,1:nDims,toplot.Null'); hold on
colorbar
ax = gca;
colormap(ax,flipud(ContextColormap))
line([tLim.start,tLim.stop],[nums.AFC,nums.AFC],'Color','black','LineWidth',1.25)
line([tLim.start,tLim.stop],[nums.AFC+nums.AW,nums.AFC+nums.AW],'Color','black','LineWidth',1.25)
line([samp,samp],[1,nDims],'Color','black','LineStyle','--','LineWidth',1.5)
xlim([tLim.start tLim.stop])
title('Null Selectivity')

nDims = size(toplot.Potent,2);
disp(num2str(nDims))
nums = PotentNums;
subplot(1,2,2)
imagesc(obj(1).time,1:nDims,toplot.Potent'); hold on
colorbar
ax = gca;
colormap(ax,flipud(ContextColormap))
line([tLim.start,tLim.stop],[nums.AFC,nums.AFC],'Color','black','LineWidth',1.25)
line([tLim.start,tLim.stop],[nums.AFC+nums.AW,nums.AFC+nums.AW],'Color','black','LineWidth',1.25)
line([samp,samp],[1,nDims],'Color','black','LineStyle','--','LineWidth',1.5)
xlim([tLim.start tLim.stop])
title('Potent Selectivity')
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
if isfield(spaceSel,'nonflip')
    toUse = spaceSel.nonflip;
else
    toUse = spaceSel;
end

% Normalize each dimension to its own max selectivity (to show dynamics of
% selectivity for each dimension over the course of the trial)
maxSelect =  max(abs(toUse),[],1);                    % Max magnitude of context selectivity for each dimension
selectNorm = [];
for c = 1:length(maxSelect)                                        % For every dimension... 
    selectNorm(:,c) = toUse(:,c)./maxSelect(c);       % Normalize all selectivity values to the max selectivity val
end 
% response = mean(selectNorm.null(start:stop,:),1);
% [~,sortix] = sort(response,'descend');
full = mean(selectNorm,1);
%[~,sortix] = sort(full,'descend');
%selectNorm = selectNorm(:,sortix);
selectNorm = selectNorm;
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

