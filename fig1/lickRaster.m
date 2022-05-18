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
% meta = loadJEB4_ALMVideo(meta,datapth); % done
% meta = loadJEB5_ALMVideo(meta,datapth); % done
% meta = loadJEB6_ALMVideo(meta,datapth); % done
meta = loadJEB7_ALMVideo(meta,datapth); % done
% meta = loadEKH1_ALMVideo(meta,datapth); % done
% meta = loadEKH3_ALMVideo(meta,datapth); % done
% meta = loadJGR2_ALMVideo(meta,datapth); % done
% meta = loadJGR3_ALMVideo(meta,datapth); % done


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

%% lick raster


% example session jeb7 4-29-2021
sessix = 2;
meta = meta(sessix);
params.probe = params.probe(sessix);
params.trialid = params.trialid{sessix};
params.cluid = params.cluid{sessix};
obj = objs{sessix};


%%
close all

clrs = getColors();

cond = [3,4];
trialOffset = 0;
for i = 1:numel(cond)
    f = figure(i); hold on;
    f.Position = [-1682         170         400         500];
    
    trix = params.trialid{cond(i)};
    trialOffset = 1;
    for j = 1:numel(trix)
%         trialOffset = trialOffset+1;
        
        
        
        goCue = obj.bp.ev.goCue(trix(j)) - obj.bp.ev.goCue(trix(j));
        sample = obj.bp.ev.sample(trix(j)) - obj.bp.ev.goCue(trix(j));
        delay = obj.bp.ev.delay(trix(j)) - obj.bp.ev.goCue(trix(j));
        
        
        lickL =  obj.bp.ev.lickL{trix(j)} - obj.bp.ev.goCue(trix(j));
        lickR =  obj.bp.ev.lickR{trix(j)} - obj.bp.ev.goCue(trix(j));
        
        
        plot([sample sample], trialOffset+[-0.5 0.5], 'k--', 'LineWidth', 2);
        plot([delay delay], trialOffset+[-0.5 0.5], 'k--', 'LineWidth', 2);
        plot([goCue goCue], trialOffset+[-0.5 0.5], 'k--', 'LineWidth', 2);
        
        if ~isempty(lickL)
            plot(lickL, trialOffset*ones(size(lickL)), '.', 'Color', clrs.lhit, 'MarkerSize',10);
            check1 = 1;
        end
        
        if ~isempty(lickR)
            plot(lickR, trialOffset*ones(size(lickR)), '.', 'Color', clrs.rhit, 'MarkerSize',10);
            check2 = 1;
        end
        
        if check1 || check2
            trialOffset = trialOffset + 1;
        end
        
        
        
    end
    xlim([-2.5 2.5])
    ylim([0 trialOffset])
    xlabel('Time (s) from go cue')
    ylabel('Trials')
    ax = gca;
    ax.YTick = [];
    ax.FontSize = 35;
    
    pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/lickraster';
    fn = ['JEB7_2021_04_29_hitmissno_notearly_' num2str(cond(i))];
    mysavefig(f,pth,fn);
    
end
















