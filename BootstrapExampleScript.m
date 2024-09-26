clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\munib\Documents\Economo-Lab\code\uninstructedMovements_v3';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
rmpath(genpath(fullfile(utilspth,'fig3/')))
rmpath(genpath(fullfile(utilspth,'fig2/')))
rmpath(genpath(fullfile(utilspth,'figx/')))
rmpath(genpath(fullfile(utilspth,'mc_stim/')))
rmpath(genpath(fullfile(utilspth,'MotionMapper/')))
rmpath(genpath(fullfile(utilspth,'musall2919/')))

% add paths for figure specific functions
addpath(genpath(pwd))

clc

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                                              % all trials       (1)
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};      % right hits, 2afc (2)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};      % left hit, 2afc   (3)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};     % right miss, 2afc (4)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};     % left miss, 2afc  (5)
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};        % 2afc hits        (6)
params.condition(end+1) = {'hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};         % aw hits          (7)
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % right hits, aw   (8)
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % left hits, aw    (9)


params.tmin = -3;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;
params.bctype = 'reflect'; % reflect, zeropad, none

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance


params.advance_movement = 0.0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

%% BOOTSTRAP

clear boot bootobj bootparams

boot.iters = 2; % number of bootstrap iterations (most papers do 1000)

% fraction of each hierarchy to sample
% (all should be 1, but can subsample by setting any to a fraction less than 1)
boot.anmfrac = 1;
boot.sessionfrac = 0.5;
boot.trialfrac = 1;

for iter = 1:boot.iters
    disp(['Iteration ' num2str(iter) '/' num2str(boot.iters)]);

    clear bootobj bootparams
    % get meta for current bootstrap iteration
    boot = hierarchicalBootstrapMeta(boot,meta,params);

    % sample data
    [bootobj,bootparams] = hierarchicalBootstrapObj(boot,obj,meta,params);

    % coding dimensions
   
    clearvars -except obj meta params boot bootobj bootparams iter allrez

    % % 2afc (early, late, go)
    cond2use = [2 3]; % left hit, right hit
    cond2proj = [2 3 4 5 8 9];
    rez_2afc = getCodingDimensions_2afc(bootobj,bootparams,cond2use,cond2proj);

    % % aw (context mode)
    cond2use = [6 7]; % hit 2afc, hit aw
    cond2proj = [2 3 4 5 8 9];
    rez_aw = getCodingDimensions_aw(bootobj,bootparams,cond2use,cond2proj);

    allrez(iter) = concatRezAcrossSessions(rez_2afc,rez_aw);

end

% CONCAT BOOTSTRAP ITERATIONS

% mean across sessions for each bootstrap iteration
temp = allrez;
for i = 1:numel(temp)
    temp(i).cd_proj = nanmean(temp(i).cd_proj,4);
    temp(i).cd_varexp = nanmean(temp(i).cd_varexp,1);
    temp(i).cd_varexp_epoch = nanmean(temp(i).cd_varexp_epoch,1);
    temp(i).selectivity_squared = nanmean(temp(i).selectivity_squared,3);
    temp(i).selexp = nanmean(temp(i).selexp,3);
end

% concatenate bootstrap iterations
rez = temp(1);
for i = 2:numel(allrez)
    rez.cd_proj = cat(4,rez.cd_proj,temp(i).cd_proj);
    rez.cd_varexp = cat(1,rez.cd_varexp,temp(i).cd_varexp);
    rez.cd_varexp_epoch = cat(1,rez.cd_varexp_epoch,temp(i).cd_varexp_epoch);
    rez.selectivity_squared = cat(3,rez.selectivity_squared,temp(i).selectivity_squared);
    rez.selexp = cat(3,rez.selexp,temp(i).selexp);
end


%% PLOTS
close all

sav = 0; % 1=save, 0=no_save

plotmiss = 0;
plotaw = 1;
plotCDProj(rez,obj(1),sav,plotmiss,plotaw,params(1).alignEvent)

% plotCDVarExp(allrez,sav)
% plotSelectivity(allrez,rez,sav)
% plotSelectivityExplained(allrez,rez,sav)

% plotCDContext_singleTrials(obj, params, rez_aw, rez_2afc)
