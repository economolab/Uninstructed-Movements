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
%% Get the reaction (after go-cue) for each trial
rt = firstJawRT(obj);
%%
cond2use = [6,11];    % All 2AFC hit trials, all AW hit trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
space = 'null_sort';
spacename = 'N_null';
[rez,me] = sortProjbyRT(rez,spacename,space,params,rt,meta,cond2use,me);

space = 'potent_sort';
spacename = 'N_potent';
[rez,me] = sortProjbyRT(rez,spacename,space,params,rt,meta,cond2use,me);
%% Group projections onto null/potent dims by RT
conds2use = [5,10];
ngroups = 4;
space = 'null_sort';
[nullgrouped] = groupProjbyRT(meta,conds2use,rez,ngroups,obj,space);
space = 'potent_sort';
[potentgrouped] = groupProjbyRT(meta,conds2use,rez,ngroups,obj,space);
%%
cols = {[0 0 0],[0.4 0.4 0.4],[0.6 0.6 0.6],[0.75 0.75 0.75]};
sm = 100;

space = 'null';                 % 'null' or 'potent'
plotRTGrouped(meta,space,nullgrouped,potentgrouped,sm,ngroups,cols)
%% Functions
function [rez,me] = sortProjbyRT(rez,spacename,space,params,rt,meta,cond2use,me)
for sessix = 1:length(meta)
    % Get the single trial projections onto all null/potent dimensions in this session
    singleproj = rez(sessix).(spacename);

    % Sort the projections for each trial according to the reaction time (separated by condition)
    for c = 1:length(cond2use)                      % For each condition...
        cond = cond2use(c);
        trix = params(sessix).trialid{cond};
        condproj = singleproj(:,trix,:);            % Get the projections on all trials of this condition
        condme = me(sessix).data(:,trix);
        condrt = rt{sessix}(trix);                  % Get the RTs for all trials of this condition

        [~,sortix] = sort(condrt,'ascend');         % Sort the RTs in ascending order
        condproj = condproj(:,sortix,:);            % Sort the projections in the same order
        condme = condme(:,sortix);
        rez(sessix).(space){c} = condproj;          % Save these sorted projections
        me(sessix).sorted{c} = condme;
    end
end
end

function [grouped] = groupProjbyRT(meta,conds2use,rez,ngroups,obj,space)
for sessix = 1:length(meta)
    for c = 1:length(conds2use)
        nTrials = size(rez(sessix).(space){c},2);
        trixPerGroup = floor(nTrials/ngroups);                  % How many trials you want to be in each group

        nDims = size(rez(sessix).(space){c},3);
        % Divide the projections onto all dimensions into nGroups
        grouped(sessix).sort{c} = NaN(length(obj(1).time),ngroups,nDims);
        cnt = 1;
        for g = 1:ngroups
            if g==ngroups
                ixrange = cnt:nTrials;
            else
                ixrange = cnt:(cnt+trixPerGroup);
            end
            grouped(sessix).sort{c}(:,g,:) = mean(rez(sessix).(space){c}(:,ixrange,:),2,'omitnan'); 
            cnt = cnt+trixPerGroup+1;
        end
    end
end
end

%%
function plotRTGrouped(meta,space,nullgrouped,potentgrouped,sm,ngroups,cols)
for sessix = 1:length(meta)
    if strcmp(space,'null')
        temp = nullgrouped;
    else
        temp = potentgrouped;
    end
    for dim = 1:size(temp(sessix).sort{1},3)
        for group = 1:ngroups
            subplot(1,2,1)
            toplot = mySmooth(temp(sessix).sort{1}(:,group,dim),sm);
            plot(toplot,'Color',cols{group},'LineWidth',1.5); hold on;
            title([space ', 2AFC'])
            ax1=gca;

            subplot(1,2,2)
            toplot = mySmooth(temp(sessix).sort{2}(:,group,dim),sm);
            plot(toplot,'Color',cols{group},'LineWidth',1.5); hold on;
            title([space ', AW'])
            ax2 = gca;
        end
        hold(ax1,'off');  hold(ax2,'off');
        sgtitle([meta(sessix).anm,' ',meta(sessix).date])
        pause
    end

end
end