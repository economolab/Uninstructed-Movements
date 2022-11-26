clear,clc,close all
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1) = {'L&hit&~stim.enable&~autowater'};             % right hits, no stim, aw off
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % left hits, no stim, aw off

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


params.advance_movement = 0.0;


%% SPECIFY DATA TO LOAD

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];

% --- ALM --- 
% meta = loadJEB6_ALMVideo(meta,datapth);
% meta = loadJEB7_ALMVideo(meta,datapth);
% meta = loadEKH1_ALMVideo(meta,datapth);
% meta = loadEKH3_ALMVideo(meta,datapth);
% meta = loadJGR2_ALMVideo(meta,datapth);
% meta = loadJGR3_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth);
% meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
meta = loadJEB14_M1TJVideo(meta,datapth);



params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA
sesh = 1;
meta = meta(sesh);
params.probe = params.probe(sesh);
[obj,params] = loadSessionData(meta,params);


%%
close all

clrs = getColors();
cols{2} = clrs.rhit;
cols{1} = clrs.lhit;


sample = mode(obj.bp.ev.sample) - mode(obj.bp.ev.(params.alignEvent));
delay = mode(obj.bp.ev.delay) - mode(obj.bp.ev.(params.alignEvent));

lw = 3;
lwx = 2;
alph = 0.5;

sm = 11;

cond = [1,2];

f1 = figure;
f1.Position = [680   730   486   248];
a1 = axes(f1);
f2 = figure;
f2.Position = [191   732   486   248];
a2 = axes(f2);

for i = 10:size(obj.psth,2)

    % psths
%     f = figure();
%     f.Position = [680   730   486   248];
    set(a1); cla(a1);
    hold(a1,'on')

    for j = 1:numel(cond)
        tempdat = mySmooth(squeeze(obj.trialdat(:,i,params.trialid{j})),sm);
        tempmean = mean(tempdat,2);
        temperr = std(tempdat,[],2) ./ sqrt(numel(params.trialid{j}));
        shadedErrorBar(obj.time,tempmean,temperr,{'Color',cols{j},'LineWidth',lw},alph, a1)

%                 temp = mySmooth(obj.psth(:,i,j),sm);
%                 plot(obj.time, temp,'Color',cols{j},'LineWidth',lw)

    end

    xline(sample,'k:','LineWidth',lwx);
    xline(delay,'k:','LineWidth',lwx);
    xline(0,'k:','LineWidth',lwx);

    xlim([obj.time(4), 2.5])
%     ylim([0 70])
    xlabel('Time (s) from go cue')
    ylabel('Firing Rate (spikes/s)')

    a1.FontSize = 14;
    hold off




    % spike raster
%     f = figure();
%     f.Position = [191   732   486   248];
    set(a2); cla(a2);
    hold(a2,'on')

    trialOffset = 1;
    clu = obj.clu{params.probe}(params.cluid(i)); 
    for j = 1:numel(cond)
        for k = 1:numel(params.trialid{j})
            trix = params.trialid{j}(k);
            ix = ismember(clu.trial,trix);

            trialtm = clu.trialtm_aligned(ix);

            scatter(trialtm,trialOffset*ones(size(trialtm)),4,'MarkerFaceColor',cols{j}, ...
            'MarkerEdgeColor','none','MarkerFaceAlpha',1)
            trialOffset = trialOffset + 1;

        end
    end

    xline(sample,'k:','LineWidth',lwx);
    xline(delay,'k:','LineWidth',lwx);
    xline(0,'k:','LineWidth',lwx);

    xlim([obj.time(4), 2.5])
    ylim([0 trialOffset])
    a2.FontSize = 14;
    hold off


    title([meta.anm ' ' meta.date ' | ' 'Cell #' num2str(i)])



    pause
end

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/example_cells';
% fn = [num2str(i) '_anmList1_sessionList1_excludeTrialTypeCount_sm_' num2str(params.smooth)];
% mysavefig(f,pth,fn);














