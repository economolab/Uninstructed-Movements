% DECODING CDlate FROM ALL KINEMATIC FEATURES
clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Hazarded Delay');
addpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\functions\HazDel funcs');

% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for                            % all trials
params.condition(1) = {'R&hit&~stim.enable&~autowater&~early'};             % R 2AFC hits, no stim
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % L 2AFC hits, no stim

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 31;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% vid features to use
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)
params.advance_movement = 0;
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
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);

%% For each session--Trials separated by delay length
plotIndiv = 'no';
for gg = 1:length(meta)
    sesh = gg;
    met = meta(sesh);

    anm = obj(sesh).pth.anm;                  % Animal name
    date = obj(sesh).pth.dt;                  % Session date
    probenum = string(met.probe);       % Which probe was used

    conditions = {1,2};

    % Find the probability of jaw [Jaw] movement at all time points in the session for trials of
    % specific conditions
    del(sesh).taxis = obj(sesh).time;
    jaw_by_cond = findJawVelocity(del(sesh).taxis, obj(sesh),conditions,'prob',params(sesh));    % (1 x conditions cell array)
    
    del(sesh).jawvel = jaw_by_cond;
end

%% Find CDs
for sesh = 1:length(meta)
    ev.sample = obj(sesh).bp.ev.sample;
    ev.delay = obj(sesh).bp.ev.delay;
    ev.goCue = obj(sesh).bp.ev.goCue;
    ev.(params(sesh).alignEvent) = obj(sesh).bp.ev.(params(sesh).alignEvent);

    % cd early mode (find using all delay lengths)
    e1 = mode(ev.sample) + 0.4 - mode(ev.(params(sesh).alignEvent));
    e2 = mode(ev.sample) + 0.8 - mode(ev.(params(sesh).alignEvent));
    times = del(sesh).taxis>e1 & del(sesh).taxis<e2;
    psth = obj(sesh).psth;
    del(sesh).cdearly_mode = calcCD_Haz(psth,times);

    % cd late mode (found only using delay length 1.2s)
    e1 = mode(ev.goCue) - 0.4 - mode(ev.(params(sesh).alignEvent));
    e2 = mode(ev.goCue) - 0.05 - mode(ev.(params(sesh).alignEvent));
    times = del(sesh).taxis>e1 & del(sesh).taxis<e2;
    psth = obj(sesh).psth;
    del(sesh).cdlate_mode = calcCD_Haz(psth,times);

    % orthogonalize
    [fns,~] = patternMatchCellArray(fieldnames(del(sesh)),{'mode'},'all');
    modes = zeros(numel(del(sesh).cdlate_mode),numel(fns));
    for i = 1:numel(fns)
        modes(:,i) = del(sesh).(fns{i});
    end

    orthModes = gschmidt(modes);

    for i = 1:numel(fns)
        del(sesh).(fns{i}) = orthModes(:,i);
    end

    % Project the PSTH for the desired delay length onto each
    % orthogonalized mode
    smooth = 51;
    for f = 1:numel(fns)
        if f==1
            tempname = strcat(fns{f}(1:7),'_latent');
        else
            tempname = strcat(fns{f}(1:6),'_latent');
        end
        for c = 1:numel(conditions)
            latent(:,c) = mySmooth(obj(sesh).psth(:,:,c)*del(sesh).(fns{f}),smooth);
        end
        del(sesh).(tempname) = latent;
    end
end
%% Avg CDs and jaw vel across sessions
clear temp
temp.cdearly_latent = NaN(length(del(1).taxis),numel(conditions),length(meta));    % time x conditions x sessions
temp.cdlate_latent = NaN(length(del(1).taxis),numel(conditions),length(meta));
temp.jawvel = NaN(length(del(1).taxis),numel(conditions),length(meta));
for sesh = 1:length(meta)
    temp.cdearly_latent(:,:,sesh) = del(sesh).cdearly_latent;
    temp.cdlate_latent(:,:,sesh) = del(sesh).cdlate_latent;
    for c = 1:numel(conditions)
        jv = mean(del(sesh).jawvel{c},2);
        temp.jawvel(:,c,sesh) = jv;
    end
end
Avg.cdearly_latent = mean(temp.cdearly_latent,3,'omitnan');  Std.cdearly = std(temp.cdearly_latent,0,3,'omitnan');
Avg.cdlate_latent = mean(temp.cdlate_latent,3,'omitnan');  Std.cdlate = std(temp.cdlate_latent,0,3,'omitnan');
Avg.jawvel = mean(temp.jawvel,3,'omitnan');  Std.jawvel = std(temp.jawvel,0,3,'omitnan');

clear temp
%%
[upperci, lowerci] = getConfInt(meta, Avg,Std,conditions);
%%
figure();
taxis = del(1).taxis;
colors = {[0 0 1],[1 0 0]};

for c = 1:numel(conditions)
     if c==1
        fn = 'R';
    else
        fn = 'L';
    end
    subplot(3,1,1)
    plot(taxis,Avg.cdearly_latent(:,c),'Color',colors{c},'LineWidth',1.5);hold on;
    patch([taxis(10:end) fliplr(taxis(10:end))],[lowerci.early.(fn)(9:end)' fliplr(upperci.early.(fn)(9:end)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
    title('CDEarly')
    ylabel('a.u.')
    xlim([-1.5 2.5])
    xline(0,'LineStyle','-.','LineWidth',1.15,'Color',[0 0 0])
    xline(0.9,'LineStyle','--','LineWidth',1,'Color',[0 0 0])
    xline(-1.3,'LineStyle','-.','LineWidth',1,'Color',[0 0 0])

    subplot(3,1,2)
    plot(taxis,Avg.cdlate_latent(:,c),'Color',colors{c},'LineWidth',1.5);hold on;
    patch([taxis(10:end) fliplr(taxis(10:end))],[lowerci.late.(fn)(9:end)' fliplr(upperci.late.(fn)(9:end)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
    title('CDLate')
    ylabel('a.u.')
    xlim([-1.5 2.5])
    xline(0,'LineStyle','-.','LineWidth',1.15,'Color',[0 0 0])
    xline(0.9,'LineStyle','--','LineWidth',1,'Color',[0 0 0])
    xline(-1.3,'LineStyle','-.','LineWidth',1,'Color',[0 0 0])

    subplot(3,1,3)
    plot(taxis,100*Avg.jawvel(:,c),'Color',colors{c},'LineWidth',1.5);hold on;
    patch([taxis(10:end) fliplr(taxis(10:end))],[100*lowerci.jawvel.(fn)(9:end)' fliplr(100*upperci.jawvel.(fn)(9:end)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
    title('Probability of jaw move')
    ylabel('%')
    xlabel(['Time from ' params(1).alignEvent ' (s)'])
    xlim([-1.5 2.5])
    xline(0,'LineStyle','-.','LineWidth',1.15,'Color',[0 0 0])
    xline(0.9,'LineStyle','--','LineWidth',1,'Color',[0 0 0])
    xline(-1.3,'LineStyle','-.','LineWidth',1,'Color',[0 0 0])
end
sgtitle('Avg across all animals')

%%
function [upperci, lowerci] = getConfInt(meta, Avg,Std,conditions)
nSessions = length(meta);
for c = 1:numel(conditions)
    if c==1
        fn = 'R';
    else
        fn = 'L';
    end
    upperci.early.(fn) = (Avg.cdearly_latent(2:end,c)+1.96*(Std.cdearly(2:end,c)/nSessions));  % Find the upper 95% confidence interval for each condition
    lowerci.early.(fn) = (Avg.cdearly_latent(2:end,c)-1.96*(Std.cdearly(2:end,c)/nSessions));  % Find lower 95% condifence interval for each condition
    %upperci.early.(fn) = fillmissing(upperci.early,'next'); lowerci.early.(fn) = fillmissing(lowerci.early,'next');

    upperci.late.(fn) = (Avg.cdlate_latent(2:end,c)+1.96*(Std.cdlate(2:end,c)/nSessions));  % Find the upper 95% confidence interval for each condition
    lowerci.late.(fn) = (Avg.cdlate_latent(2:end,c)-1.96*(Std.cdlate(2:end,c)/nSessions));  % Find lower 95% condifence interval for each condition
    %upperci.late.(fn) = fillmissing(upperci.late,'next'); lowerci.late.(fn) = fillmissing(lowerci.late,'next');

    upperci.jawvel.(fn) = (Avg.jawvel(2:end,c)+1.96*(Std.jawvel(2:end,c)/nSessions));  % Find the upper 95% confidence interval for each condition
    lowerci.jawvel.(fn) = (Avg.jawvel(2:end,c)-1.96*(Std.jawvel(2:end,c)/nSessions));  % Find lower 95% condifence interval for each condition
    %upperci.jawvel.(fn) = fillmissing(upperci.jawvel,'next'); lowerci.jawvel.(fn) = fillmissing(lowerci.jawvel,'next');
end
end