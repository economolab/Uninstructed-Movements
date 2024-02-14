% Finding CDContext from neural activity that resides within the Null and Potent spaces
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

load([basepth '\Uninstructed-Movements\ContextColormap.mat']);
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
    [blockid,nBlocks] = getBlockNum_AltContextTask(sessix,obj);
    
    % -- input data
     trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix));
    zscored(sessix).trialdat =  trialdat_zscored;

    % -- Calculate the null and potent spaces for each session
    cond2use = [2 3 4 5];   % All 2AFC hit trials, all AW hit trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
    nullalltime = 0;        % use all time points to estimate null space if 1
    AWonly = 0;             % use only AW to find null and potent spaces 
    delayOnly = 0;          % use only delay period to find null and potent spaces
    cond2proj = [2 3];       % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj,nullalltime,AWonly,delayOnly);

    % -- Find coding dimensions from RECONSTRUCTED full neural activity which is reconstructed from the null and potent spaces
    cond2use = [1 2];            % (NUMBERING ACCORDING TO THE CONDITIONS PROJECTED INTO NULL AND POTENT SPACES, i.e. which of the conditions specified in 'cond2proj' above do you want to use?)
    cond2proj = [1 2];           % 2AFC hits, AW hits, 2AFC miss, AW miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3];   % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    cd_null(sessix) = getCodingDimensions_Context_NonStation(rez(sessix).recon_psth.null,...
        rez(sessix).recon.null,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,nBlocks,blockid);
    
    cd_potent(sessix) = getCodingDimensions_Context_NonStation(rez(sessix).recon_psth.potent,...
        rez(sessix).recon.potent,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj,nBlocks,blockid);

end
%% Project single trials onto Null and Potent CDs
disp('----Projecting single trials onto CDContext----')
cd = 'context';

[cd_null,cd_potent] = getNPSingleTrialProjs(obj,cd,cd_null,cd_potent,rez); 
%% Plot heatmap of single trialCDContext Null and Potent over the course of a session
% tmax = 0;                               % what you want max of xlim to be
% for sessix = 1:length(meta)%sessrange
%     plotNP_CDCont_Heatmap(sessix, cd_null, cd_potent,obj,tmax)
% end
% 
% %% Find trials in which the animal switches between contexts
% for sessix = 1:numel(meta)
%     bp = obj(sessix).bp;
%     [toAW_ix, toAFC_ix] = findSwitchTrials(bp);
%     % toAW_ix = the first trial in an AW block; toAFC_ix = the first trial in a 2AFC block
%     obj(sessix).toAW_ix = toAW_ix;  obj(sessix).toAFC_ix = toAFC_ix;
% end
% %% Find avg context mode across trials in a given session on switch trials
% clear temp
% nBufferTrials = 10;                              % Number of trials that you want to look at before and after switches
% % Time period over which you want to average CDContext
% trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
% start = find(obj(1).time>trialstart,1,'first');
% samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
% stop = find(obj(1).time<samp,1,'last');
% 
% for s=1:2
%     if s==1
%         space = cd_null;
%         spacename = 'Null';
%     else
%         space = cd_potent;
%         spacename = 'Potent';
%     end
%     for sessix = 1:numel(meta)
%         CDCont = space(sessix).singleProj.context;
% 
%         % Find avg CDContext on 2AFC --> AWswitch trials (during the presample period)
%         switchtype = 'toAW_ix';
%         contswitch(sessix).toAW_CDCont.(spacename) = findCDCont_SwitchAligned(nBufferTrials, obj, sessix, CDCont,switchtype,start,stop);
% 
%         % Find avg CDContext on AW --> 2AFC switch trials (during the presample period)
%         switchtype = 'toAFC_ix';
%         contswitch(sessix).to2AFC_CDCont.(spacename) = findCDCont_SwitchAligned(nBufferTrials, obj, sessix, CDCont,switchtype,start,stop);
%     end
% end
% %%
% % Normalize the CDContext values to the max of the absolute values (for
% % each session--i.e. each session will have CDCont values between -1 and 1)
% % for sessix = 1:length(meta)
% %     blah = contswitch(sessix);
% %     blah = normalizeCDCont(blah);
% %     contswitch(sessix) = blah;
% % end
% %% Plot switch-aligned CDContext
% % Concatenate all switch-aligned CDContexts from across sessions
% for s=1:2
%     if s==1
%         spacename = 'Null';
%     else
%         spacename = 'Potent';
%     end
%     toAW = []; to2AFC = [];
%     for sessix = 1:length(meta)
%         toAW = [toAW; contswitch(sessix).toAW_CDCont.(spacename)];
%         to2AFC = [to2AFC; contswitch(sessix).to2AFC_CDCont.(spacename)];
%     end
%     col = [0.35 0.35 0.35];
%     nSessions = length(meta);
%     tRange = -nBufferTrials:nBufferTrials;           % Number of trials that you want to plot
%     stdCD = getStd(toAW,to2AFC);                     % Get the standard deviation of CDCont across all trials (from all sessions)
%     alpha = 0.2;                                     % Opacity of confidence intervals
% 
%     plotSwitchAlignedCDCont(toAW, to2AFC,stdCD, nSessions, tRange, alpha, col,nBufferTrials)
%     sgtitle(spacename)
% end
%% Concatenate coding dimensions in the null and potent space across sessions
cd_null_all = concatRezAcrossSessions_Context(cd_null);
cd_potent_all = concatRezAcrossSessions_Context(cd_potent);
%% plots

close all

sav = 0;

% -----------------------------------------------------------------------
% -- Null and Potent Space Single Trial Projections --
% -----------------------------------------------------------------------

% % - projections showing move / quiet somehow
cond2use = [2 3];                       % 2AFC hits, AW hits (reference to params.condition)
ndims = 4;                              % top ndims variance explaining dimensions
%plotSingleTrialNPHeatmaps(rez,params,me,ndims,cond2use,meta);

%%
% -----------------------------------------------------------------------
% -- Null and Potent Space Trial-averaged Projections --
% -----------------------------------------------------------------------
ndims = 5;                              % how many n/p dimensions to plot, in order of variance explained
cond2plot = 1:2;                        % 2AFC hit, AW hit (reference to conditions used to find null_psth)
% plot_NP_PSTH(rez,obj,params,ndims,cond2plot,meta)
% plotTopDims_NP_Proj(meta,rez,obj,cond2plot,ndims)
%%
% -----------------------------------------------------------------------
% -- Coding Dimensions --
% -----------------------------------------------------------------------
plotmiss = 0;

figure();

titlestring = 'Null';
subplot(1,2,1)
plotCDProj_Context(cd_null_all,cd_null,sav,titlestring,plotmiss)

titlestring = 'Potent';
subplot(1,2,2)
plotCDProj_Context(cd_potent_all,cd_potent,sav,titlestring,plotmiss)

%plotCDContext_Selectivity(cd_null_all, cd_potent_all,cd_null)

% plotCDContext_SelectivityScatter(cd_null_all, cd_potent_all,cd_potent,obj(1),params(1))


%plotNP_CD_Context(cd_null_all,cd_null,cd_potent_all,cd_potent)
%%

function stdCD = getStd(toAW,to2AFC)

stdCD.toAW = std(toAW,0,1,'omitnan');     % Get standard deviation of switch aligned presamp CDContext across all trials
stdCD.to2AFC = std(to2AFC,0,1,'omitnan');
end

function blah = normalizeCDCont(blah)
maxblah = max(abs(blah.toAW_CDCont)); blah.toAW_CDCont = blah.toAW_CDCont./maxblah;
maxblah = max(abs(blah.to2AFC_CDCont)); blah.to2AFC_CDCont = blah.to2AFC_CDCont./maxblah;
end

% Plotting
function plotSwitchAlignedCDCont(toAW, to2AFC,stdCD, nSessions, tRange, alpha, col,nBufferTrials)
figure();
subplot(1,2,1)
toplot = mean(toAW,1,'omitnan');
plot(tRange,toplot)
hold on;
ax = gca;
err = 1.96*(stdCD.toAW/nSessions);
shadedErrorBar(tRange, toplot, err ,{'Color',col,'LineWidth',2}, alpha, ax)
xline(0,'k--')
title('2AFC to AW switches')
xlabel('Trials to context switch')
xlim([-nBufferTrials, nBufferTrials]);
ylabel('Normalized MOVE-CDContext proj (a.u.)')

subplot(1,2,2)
toplot = mean(to2AFC,1,'omitnan');
plot(tRange,toplot)
hold on;
ax = gca;
err = 1.96*(stdCD.toAW/nSessions);
shadedErrorBar(tRange, toplot, err ,{'Color',col,'LineWidth',2}, alpha, ax)
xline(0,'k--')
title('AW to 2AFC switches')
xlabel('Trials to context switch')
xlim([-nBufferTrials, nBufferTrials]);
ylabel('Normalized MOVE-CDContext proj (a.u.)')
end

