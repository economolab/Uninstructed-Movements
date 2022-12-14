clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig1/')))
rmpath(genpath(fullfile(utilspth,'mc_stim/')))

% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % error left, no stim, aw off

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
                        {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0.0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];
% 
% meta = loadJEB6_ALMVideo(meta,datapth);
% meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth);
% meta = loadEKH3_ALMVideo(meta,datapth);
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);

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
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat);

    % -- null and potent spaces
    cond2use = [2 3 4 5]; % right hit, left hit, right miss, left miss
    first = 'null'; % 'null'  'potent'
    rez(sessix) = singleTrial_pca_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, first);

    % -- coding dimensions
    cond2use = [1 2]; % right hits, left hits (corresponding to null/potent psths in rez)
    cond2proj = [1:4]; % right hits, left hits, right miss, left miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
    cd_null(sessix) = getCodingDimensions(rez(sessix).N_null_psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat,cond2proj);
    cd_potent(sessix) = getCodingDimensions(rez(sessix).N_potent_psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat,cond2proj);
end


% concatenate coding dimension results across sessions
cd_null_all = concatRezAcrossSessions(cd_null);
cd_potent_all = concatRezAcrossSessions(cd_potent);


%% plots

close all

sav = 0;


% -----------------------------------------------------------------------
% -- Null and Potent Space Single Trial Projections --
% -----------------------------------------------------------------------
 
% % - projections showing move / quiet somehow
cond2use = [2 3]; % right hits, left hits
ndims = 4; % top ndims variance explaining dimensions
% plotSingleTrialNPHeatmaps(rez,params,mendims,cond2use,meta);


% % - how much variance in move and non-move time points
cond2use = [2 3]; % right hits, left hits
% plotVarianceInEpochs(rez,me,params,cond2use);                       

% % - ve
% plotVarianceExplained_NP(rez);

% % - ve over time (TODO)
% % % plotVarianceExplained_NP_overTime(rez);


% -----------------------------------------------------------------------
% -- Null and Potent Space Trial-averaged Projections --
% -----------------------------------------------------------------------

ndims = 5; % how many n/p dimensions to plot, in order of variance explained
cond2plot = [1 2]; % right hit, left hit
% plot_NP_PSTH(rez,obj,params,ndims,cond2plot,meta)


% -----------------------------------------------------------------------
% -- Coding Dimensions --
% -----------------------------------------------------------------------

titlestring = 'Null';
% plotCDProj(cd_null_all,cd_null,sav,titlestring)
% plotCDVarExp(cd_null_all,sav,titlestring)
% plotSelectivity(cd_null_all,cd_null,sav,titlestring)
% plotSelectivityExplained(cd_null_all,cd_null,sav,titlestring)

titlestring = 'Potent';
% plotCDProj(cd_potent_all,cd_potent,sav,titlestring)
% plotCDVarExp(cd_potent_all,sav,titlestring)
% plotSelectivity(cd_potent_all,cd_potent,sav,titlestring)
% plotSelectivityExplained(cd_potent_all,cd_potent,sav,titlestring)


%% t=0 is the go cue, but only on trials where the animals were not moving PRIOR to the go cue
% same plots as plotSelectivityExplained

%% t=0 is transitions between non-movement and movement that do not coincide with the go cue
% same plots as plotSelectivityExplained

stitch_dist = 0.05; % in seconds, stitch together movement bouts shorter than this % 0.025
purge_dist = 0.1; % in seconds, remove move bouts shorter than this value, after stitching complete % 0.1
tbout = 0.3; % move/non-move bout/transition required to be at least this long in seconds % 0.3
[dat.mdat,dat.mdat_leg,dat.qdat,dat.qdat_leg,newme] = nonGC_moveTransitions(obj,me,params,stitch_dist,purge_dist,tbout);

%%
close all

sav = 0;

% TODO: rewrite these functions to be more general

dim = 1; % dim to plot (most to least variance explained)

% plot_nonGC_moveTransitions_singleTrials(dat,obj,newme,rez,params,dim,meta)

% plot_nonGC_moveTransitions_singleTrials_v2(dat,obj,newme,rez,params,dim,meta) % random sampling of move to quiet bouts
% plot_nonGC_moveTransitions_singleTrials_v3(dat,obj,newme,rez,params,dim,meta) % random sampling of quiet to move bouts

ndims = 10;
% plot_nonGC_moveTransitions_singleTrials_v4(dat,obj,newme,rez,params,ndims,meta,sav) % all bouts heat map, move to quiet 
% plot_nonGC_moveTransitions_singleTrials_v5(dat,obj,newme,rez,params,ndims,meta,sav) % all bouts heat map, quiet to move

ndims = 10;
% plot_nonGC_moveTransitions_trialAvg(dat,obj,newme,rez,params,meta,ndims) % separate tiles for each dimension

% plot_nonGC_moveTransitions_trialAvg_v2(dat,obj,newme,rez,params,meta,ndims) % all dimensions plotted on same axis, 

% plot_nonGC_moveTransitions_trialAvg_v3(dat,obj,newme,rez,params,meta,ndims) % mean,stderr across dimensions

% plot_nonGC_moveTransitions_trialAvg_v4(dat,obj,newme,rez,params,meta,ndims) % var across dimensions, move to quiet
% plot_nonGC_moveTransitions_trialAvg_v5(dat,obj,newme,rez,params,meta,ndims) % var across dimensions, quiet to move

plot_nonGC_moveTransitions_trialAvg_v6(dat,obj,newme,rez,params,meta,ndims) % sumsqmag across dimensions, move to quiet
plot_nonGC_moveTransitions_trialAvg_v7(dat,obj,newme,rez,params,meta,ndims) % sumsqmag across dimensions, quiet to move



%% cross corr me and potent space
clear dat

maxlag = int64( 5 / params(1).dt );


for sessix = 1:numel(rez)
    dat.null.ts = rez(sessix).N_null;
    dat.potent.ts = rez(sessix).N_potent;
    dat.me = zscore(me(sessix).data);
    
    for dimix = 1:size(dat.null.ts,3)
        for trix = 1:size(dat.null.ts,2)
            [dat.null.xc{sessix}(:,trix,dimix),dat.lags] = xcorr(dat.me(:,trix), dat.null.ts(:,trix,dimix), maxlag, 'normalized');
            [dat.potent.xc{sessix}(:,trix,dimix),dat.lags] = xcorr(dat.me(:,trix), dat.potent.ts(:,trix,dimix), maxlag, 'normalized');
        end
    end
    
end


dat.null.xc_mean_across_trials = cellfun(@(x) squeeze(nanmean(x,2)), dat.null.xc ,'UniformOutput',false);
dat.potent.xc_mean_across_trials = cellfun(@(x) squeeze(nanmean(x,2)), dat.potent.xc ,'UniformOutput',false);

dat.null.xc_mean_across_dims = cellfun(@(x) squeeze(nanmean(x,2)),dat.null.xc_mean_across_trials,'UniformOutput',false);
dat.potent.xc_mean_across_dims = cellfun(@(x) squeeze(nanmean(x,2)),dat.potent.xc_mean_across_trials,'UniformOutput',false);

dat.null.xc_meanall = cell2mat(dat.null.xc_mean_across_dims);
dat.potent.xc_meanall = cell2mat(dat.potent.xc_mean_across_dims);

%% plot lags against neurons for trial-averaged data
close all


tm = dat.lags .* params(1).dt;
ix1 = find(tm>=-1,1,'first');
ix2 = find(tm<=1,1,'last');

% null

f = figure;
t = tiledlayout('flow');
for i = 1:numel(dat.null.xc_mean_across_trials)
    ax = nexttile;
    temp = dat.null.xc_mean_across_trials{i};
    [~,ix] = sort(mean(temp(ix1:ix2,:),1));
    
    imagesc(tm ,1:size(temp,2),temp(:,ix)')
    title([meta(i).anm '_' meta(i).date],'Interpreter','none')
end
xlabel(t,'lags (s)')
ylabel(t,'dims')
title(t,'Null, cross correlation b/w motion energy and single trial fr')

% potent

f = figure;
t = tiledlayout('flow');
for i = 1:numel(dat.potent.xc_mean_across_trials)
    ax = nexttile;
    temp = dat.potent.xc_mean_across_trials{i};
    [~,ix] = sort(mean(temp(ix1:ix2,:),1));
    
    imagesc(tm ,1:size(temp,2),temp(:,ix)')
    title([meta(i).anm '_' meta(i).date],'Interpreter','none')
end
xlabel(t,'lags (s)')
ylabel(t,'dims')
title(t,'Potent, cross correlation b/w motion energy and single trial fr')

%% mean lag for all dims

close all

xlims = [-1 1];

% null

tm = dat.lags .* params(1).dt;

rr_mean = mean(dat.null.xc_meanall,2,'omitnan');
rr_err = std(dat.null.xc_meanall,[],2,'omitnan') ./ numel(meta);

[maxs,imaxs] = max(abs(dat.null.xc_meanall)); 
tmaxs = tm(imaxs);

meanalph = 0.5;
linealph = 0.8;
lw = 2;
cols = linspecer(size(dat.null.xc_meanall,2));

f = figure; hold on;
ax = gca;
% ax = subplot(2,1,1);
% hold on;
% for i = 1:size(dat.null.xc_meanall,2)
%     patchline(tm,dat.null.xc_meanall(:,i),'EdgeColor',cols(i,:),'LineWidth',1,'EdgeAlpha',linealph)
%     scatter(tmaxs(i),dat.null.xc_meanall(imaxs(i),i),60,'MarkerEdgeColor','none','MarkerFaceColor',cols(i,:))
% end
shadedErrorBar(tm,rr_mean,rr_err,{'Color',[62, 168, 105]./255,'LineWidth',lw},meanalph,ax)
% xlim(xlims)
xlabel('Time shift (s)')
ylabel('correlation')

% title('null')

% potent

tm = dat.lags .* params(1).dt;

rr_mean = mean(dat.potent.xc_meanall,2,'omitnan');
rr_err = std(dat.potent.xc_meanall,[],2,'omitnan') ./ numel(meta);

[maxs,imaxs] = max(abs(dat.potent.xc_meanall)); 
tmaxs = tm(imaxs);

meanalph = 0.5;
linealph = 0.8;
lw = 2;
cols = linspecer(size(dat.potent.xc_meanall,2));

% f = figure; 
% ax = gca;
% ax = subplot(2,1,1);
% hold on;
% for i = 1:size(dat.potent.xc_meanall,2)
%     patchline(tm,dat.potent.xc_meanall(:,i),'EdgeColor',cols(i,:),'LineWidth',1,'EdgeAlpha',linealph)
%     scatter(tmaxs(i),dat.potent.xc_meanall(imaxs(i),i),60,'MarkerEdgeColor','none','MarkerFaceColor',cols(i,:))
% end
shadedErrorBar(tm,rr_mean,rr_err,{'Color',[255, 56, 140]./255,'LineWidth',lw},meanalph,ax)
% xlim(xlims)
xlabel('Time shift (s)')
ylabel('correlation')


% title('potent')

% 
% % potent - null
% 
% tm = dat.lags .* params(1).dt;
% 
% pn_minus = dat.potent.xc_meanall - dat.null.xc_meanall;
% 
% rr_mean = mean(pn_minus,2,'omitnan');
% rr_err = std(pn_minus,[],2,'omitnan') ./ numel(meta);
% 
% [maxs,imaxs] = max(abs(pn_minus)); 
% tmaxs = tm(imaxs);
% 
% meanalph = 0.5;
% linealph = 0.2;
% lw = 2;
% cols = linspecer(size(pn_minus,2));
% 
% f = figure; 
% ax = gca;
% shadedErrorBar(tm,rr_mean,rr_err,{'Color','k','LineWidth',lw},meanalph,ax)
% % xlim([-0.1 0.1])
% xlabel('Time shift (s)')
% ylabel('correlation')
% 
% title('potent - null')


%% cross corr me and potent space (ssm)
clear dat

maxlag = int64( 5 / params(1).dt );


for sessix = 1:numel(rez)
    dat.null.ts = sum(rez(sessix).N_null.^2,3);
    dat.potent.ts = rez(sessix).N_potent;
    dat.me = zscore(me(sessix).data);

    for trix = 1:size(dat.null.ts,2)
        [dat.null.xc{sessix}(:,trix),dat.lags] = xcorr(dat.me(:,trix), dat.null.ts(:,trix), maxlag, 'normalized');
        [dat.potent.xc{sessix}(:,trix),dat.lags] = xcorr(dat.me(:,trix), dat.potent.ts(:,trix), maxlag, 'normalized');
    end

end


dat.null.xc_mean_across_trials = cellfun(@(x) squeeze(nanmean(x,2)), dat.null.xc ,'UniformOutput',false);
dat.potent.xc_mean_across_trials = cellfun(@(x) squeeze(nanmean(x,2)), dat.potent.xc ,'UniformOutput',false);


dat.null.xc_meanall = cell2mat(dat.null.xc_mean_across_trials);
dat.potent.xc_meanall = cell2mat(dat.potent.xc_mean_across_trials);


%%
close all

xlims = [-1 1];

% null

tm = dat.lags .* params(1).dt;

rr_mean = mean(dat.null.xc_meanall,2,'omitnan');
rr_err = std(dat.null.xc_meanall,[],2,'omitnan') ./ numel(meta);

[maxs,imaxs] = max(abs(dat.null.xc_meanall)); 
tmaxs = tm(imaxs);

meanalph = 0.5;
linealph = 0.8;
lw = 2;
cols = linspecer(size(dat.null.xc_meanall,2));

f = figure; hold on;
ax = gca;
% ax = subplot(2,1,1);
% hold on;
% for i = 1:size(dat.null.xc_meanall,2)
%     patchline(tm,dat.null.xc_meanall(:,i),'EdgeColor',cols(i,:),'LineWidth',1,'EdgeAlpha',linealph)
%     scatter(tmaxs(i),dat.null.xc_meanall(imaxs(i),i),60,'MarkerEdgeColor','none','MarkerFaceColor',cols(i,:))
% end
shadedErrorBar(tm,rr_mean,rr_err,{'Color',[62, 168, 105]./255,'LineWidth',lw},meanalph,ax)
xlim(xlims)
xlabel('Time shift (s)')
ylabel('correlation')

% title('null')

% potent

tm = dat.lags .* params(1).dt;

rr_mean = mean(dat.potent.xc_meanall,2,'omitnan');
rr_err = std(dat.potent.xc_meanall,[],2,'omitnan') ./ numel(meta);

[maxs,imaxs] = max(abs(dat.potent.xc_meanall)); 
tmaxs = tm(imaxs);

meanalph = 0.5;
linealph = 0.8;
lw = 2;
cols = linspecer(size(dat.potent.xc_meanall,2));

% f = figure; 
% ax = gca;
% ax = subplot(2,1,1);
% hold on;
% for i = 1:size(dat.potent.xc_meanall,2)
%     patchline(tm,dat.potent.xc_meanall(:,i),'EdgeColor',cols(i,:),'LineWidth',1,'EdgeAlpha',linealph)
%     scatter(tmaxs(i),dat.potent.xc_meanall(imaxs(i),i),60,'MarkerEdgeColor','none','MarkerFaceColor',cols(i,:))
% end
shadedErrorBar(tm,rr_mean,rr_err,{'Color',[255, 56, 140]./255,'LineWidth',lw},meanalph,ax)
xlim(xlims)
xlabel('Time shift (s)')
ylabel('correlation')










