clear,clc,close all

addpath(genpath(pwd))

% finds cd early, late, go as defined in economo 2018


%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
params.condition(end+1) = {'(hit|miss|no)&R&~autowater'}; % right trials
params.condition(end+1) = {'(hit|miss|no)&L&~autowater'}; % left trials


% set conditions used for finding activity modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %rhit, aw off 
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};  %lhit, aw off 
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %rmiss, aw off 
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']}; %lmiss, aw off 
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};    % hit, aw off 

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'Fair','Good','Great','Excellent','single'}; 

%% SET METADATA

datapth = '/Users/Munib/Documents/Economo-Lab/data/';

meta = [];
meta = loadJEB4_ALMVideo(meta,datapth); % done
meta = loadJEB5_ALMVideo(meta,datapth); % done
meta = loadJEB6_ALMVideo(meta,datapth); % done
meta = loadJEB7_ALMVideo(meta,datapth); % done
meta = loadEKH1_ALMVideo(meta,datapth); % done
meta = loadEKH3_ALMVideo(meta,datapth); % done
meta = loadJGR2_ALMVideo(meta,datapth); % done
meta = loadJGR3_ALMVideo(meta,datapth); % done


params.probe = [meta.probe]; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

%% LOAD AND PROCESS DATA

objs = loadObjs(meta);


for metaix = 1:numel(meta)
    obj = objs{metaix};
    disp('______________________________________________________')
    disp(['Processing data for session ' [meta(metaix).anm '_' meta(metaix).date]])
    disp(' ')
    [sessparams{metaix},sessobj{metaix}] = processData(obj,params,params.probe(metaix));
end

% clean up sessparams and sessobj
for metaix = 1:numel(meta)
    params.trialid{metaix} = sessparams{metaix}.trialid;
    params.cluid{metaix} = sessparams{metaix}.cluid{params.probe(metaix)};
    
    objs{metaix} = sessobj{metaix};
    objs{metaix}.psth = objs{metaix}.psth{params.probe(metaix)};
    objs{metaix}.trialdat = objs{metaix}.trialdat{params.probe(metaix)};
    objs{metaix}.presampleFR = objs{metaix}.presampleFR{params.probe(metaix)};
    objs{metaix}.presampleSigma = objs{metaix}.presampleSigma{params.probe(metaix)};
end

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')

%% COMBINE DATA ACROSS SESSIONS

% TODO: Don't normalize here, and see how this changes cds and selectivity.
% currently, total selectivity looks weird with normalized data (no
% ramping)

% psth = normalizePSTH(objs{1});
psth = objs{1}.psth;
ev.sample = objs{1}.bp.ev.sample;
ev.delay = objs{1}.bp.ev.delay;
ev.goCue = objs{1}.bp.ev.goCue;
ev.(params.alignEvent) = objs{1}.bp.ev.(params.alignEvent);
for sessix = 2:numel(objs)
    obj = objs{sessix};
    
%     temppsth = normalizePSTH(obj);
    temppsth = obj.psth;
    
    psth = cat(2,psth,temppsth);
    
    ev.sample = [ev.sample ; obj.bp.ev.sample];
    ev.delay = [ev.delay ; obj.bp.ev.delay];
    ev.goCue = [ev.goCue ; obj.bp.ev.goCue];
    ev.(params.alignEvent) = [ev.(params.alignEvent) ; obj.bp.ev.(params.alignEvent)];
    
end

psth(isnan(psth)) = 0;
psth(isinf(psth)) = 0;


psth = psth(:,:,[1,2]);

%%
close all

clrs = getColors();
cols{1} = clrs.rhit;
cols{2} = clrs.lhit;


sample = mode(ev.sample) - mode(ev.(params.alignEvent));
delay = mode(ev.delay) - mode(ev.(params.alignEvent));

lw = 6;
lwx = 2;

sm = 101;

cond = [1,2];

f = figure();
f.Position = [-1351         382         483         353];
for i = 1:size(psth,2)
    clf(f); hold on
    for j = 1:numel(cond)
        
        temp = mySmooth(psth(:,i,j),sm);
        plot(objs{1}.time, temp,'Color',cols{j},'LineWidth',lw)
        
    end
    
    xline(sample,'k--','LineWidth',lwx);
    xline(delay,'k--','LineWidth',lwx); 
    xline(0,'k--','LineWidth',lwx); 
    
    xlim([-2.45 2.5])
    xlabel('Time (s) from go cue')
    ylabel('Firing Rate (spikes/s)')
    title(['Cell ' num2str(i)])
    ax = gca;
    ax.FontSize = 30;
    hold off
    pause
end

% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/example_cells';
% fn = [num2str(i) '_anmList1_sessionList1_excludeTrialTypeCount_sm_' num2str(params.smooth)];
% mysavefig(f,pth,fn);














