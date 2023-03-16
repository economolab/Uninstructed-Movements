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
    cd_context(sessix) =getCodingDimensions_Context(obj(sessix).psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);
end
%% Sanity check
% for sessix = 1:length(meta)
%     for ii = 1:3
%         switch ii
%             case 1
%                 space = cd_context(sessix);
%                 titl = 'Full';
%             case 2
%                 space = cd_null(sessix);
%                 titl = 'Null';
%             case 3
%                 space = cd_potent(sessix);
%                 titl = 'Potent';
%         end
%         subplot(3,1,ii)
%         for c = 1:2
%             plot(space.cd_proj(:,c)); hold on;
%         end
%         hold off;
%         title(titl)
%     end
%     pause
% end
%%
clearvars -except obj meta params me sav rez cd_null cd_potent cd_context
% Find the times corresponding to trial start and the sample period
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>trialstart,1,'first');
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');
%% Find correlation/selectivity explained by CDContext found from the reconstructed PSTHs from null or potent spaces
% Asking about the selectivty explained by the null/potent spaces only in the presample period 
sel_explained.null = NaN(1,length(meta));
sel_explained.potent = NaN(1,length(meta));
for sessix = 1:length(meta)
    fullsel = cd_context(sessix).cd_proj(:,1)-cd_context(sessix).cd_proj(:,2);
    nullsel = cd_null(sessix).cd_proj(:,1)-cd_null(sessix).cd_proj(:,2);
    potentsel = cd_potent(sessix).cd_proj(:,1)-cd_potent(sessix).cd_proj(:,2);

    nullR2 = corrcoef(fullsel(start:stop),nullsel(start:stop));
    nullR2 = nullR2(2);

    potentR2 = corrcoef(fullsel(start:stop),potentsel(start:stop));
    potentR2 = potentR2(2);

    sel_explained.null(sessix) = abs(nullR2);
    sel_explained.potent(sessix) = abs(potentR2);
end
%% Bar plot of R2 values on each session between true and predicted Null and Potent CDlate
figure();
colors = getColors_Updated();
X = categorical({'Null','Potent'});
X = reordercats(X,{'Null','Potent'});
Y = [mean(sel_explained.null), mean(sel_explained.potent)];
b = bar(X,Y); hold on;
b.FaceColor = 'flat';
b.FaceAlpha = 0.5;
b.CData(1,:) = colors.null;
b.CData(2,:) = colors.potent;
scatter(1,sel_explained.null,35,[0 0 0],'filled','MarkerEdgeColor','black')
scatter(2,sel_explained.potent,35,[0 0 0],'filled','MarkerEdgeColor','black')
for sessix = 1:length(sel_explained.null)
    plot([1,2],[sel_explained.null(sessix),sel_explained.potent(sessix)],'Color','black')
end

sigcutoff = 0.05;
hyp = ttest(sel_explained.null,sel_explained.potent,'Alpha',sigcutoff);
if hyp
    scatter(1.5,1.2,30,'*','MarkerEdgeColor','black')
    title(['Presamp selectivity explained by null/pot CDCont; p<' num2str(sigcutoff)])
else
    title(['Presamp selectivity explained by null/pot CDCont; p>' num2str(sigcutoff)])
end
ylabel("R^2 value")

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
[~,sortix] = sort(full,'descend');
selectNorm = selectNorm(:,sortix);
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

