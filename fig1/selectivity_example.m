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
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};          % right hits, no stim, aw on
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};          % left hits, no stim, aw on
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};        % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};        % error left, no stim, aw off
params.condition(end+1) = {'~hit&~miss&~stim.enable&~autowater&~early'};    % ignore
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};           % hit 2afc
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};            % hit aw


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
params.smooth = 51;

% cluster qualities to use
% params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality
params.quality = {'all'}; 

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

%% Remove unwanted sessions

% remove sessions if:
% 1) less than 40 trials of rhit and lhit each (same as
% 2) atleast 10 clusters of quality that is not noise
use = false(size(objs));
for i = 1:numel(use)
    check(1) = numel(params.trialid{i}{1}) > 40;
    check(2) = numel(params.trialid{i}{2}) > 40;
    check(3) = numel(params.cluid{i}) >= 10;
    if all(check)
        use(i) = true;
    end
end

meta = meta(use);
objs = objs(use);
params.probe = params.probe(use);
params.trialid = params.trialid(use);
params.cluid = params.cluid(use);


% % just one session
% use = 2;
% meta = meta(use);
% obj = objs(use);
% params.probe = params.probe(use);
% params.trialid = params.trialid(use);
% params.cluid = params.cluid(use);
% 
% obj = objs{1};

%%

close all

clrs = getColors();
lw = 4;
alph = 0.5;
sm = 31;

f = figure;
for i = 1:numel(params.cluid{1})
    clf
    clu = params.cluid{1}(1);
    right = mySmooth(squeeze(obj.psth(:,i,1)),sm);
    left = mySmooth(squeeze(obj.psth(:,i,2)),sm);
    
    
    hold on
    plot(obj.time,right - left, 'Color', [224, 78, 237]./255, 'LineWidth',lw)
    patchline(obj.time,right,'EdgeColor',clrs.rhit,'EdgeAlpha',alph,'LineWidth',2)
    patchline(obj.time,left,'EdgeColor',clrs.lhit,'EdgeAlpha',alph,'LineWidth',2)
    yline(0,'k--','LineWidth',2)
    xlabel('Time (s) from go cue')
    ylabel('Spikes / sec')
    xlim([obj.time(15) obj.time(end)]);
    ax = gca;
    ax.FontSize = 20;
    
    
    
    
    pause
    
end

%%

sm = 21;

selpsth = objs{1}.psth;
for sessix = 2:numel(objs)
    obj = objs{sessix};    
    selpsth = cat(2,selpsth,obj.psth);
end
selpsth(isnan(selpsth)) = 0;
selpsth(isinf(selpsth)) = 0;

sel = mySmooth(selpsth(:,:,1) - selpsth(:,:,2),sm);

% sel(sel>abs(min(min(sel)))) = abs(min(min(sel)));

%%

close all

mask = obj.time > -1.4 & obj.time < -1.0;
means = mean(sel(mask,:),1);
[~,ix] = sort(means);
selsorted = sel(:,ix);

f = figure;
imagesc(obj.time,1:size(selsorted,2),selsorted')

% cmap = flipud(colormap(hot));
cmap = colormap(hot);

colormap(cmap)
c = colorbar;
caxis([-40 40])
ax = gca;
ax.FontSize = 20;
c.Label.String = 'Selectivity (spks/sec)';
xlabel('Time (s) from go cue')
ylabel('Neurons')

%%
close all
mask = obj.time > -0.5 & obj.time < -0.1;
means = mean(sel(mask,:),1);
[~,ix] = sort(means);
selsorted = sel(:,ix);

f = figure;
imagesc(obj.time,1:size(selsorted,2),selsorted')

% cmap = flipud(colormap(hot));
cmap = colormap(hot);

colormap(cmap)
c = colorbar;
caxis([-40 40])
ax = gca;
ax.FontSize = 20;
c.Label.String = 'Selectivity (spks/sec)';
xlabel('Time (s) from go cue')
ylabel('Neurons')


f = figure; 
imagesc(obj.time(400),1:size(selsorted,2),selsorted(400,:)');
cmap = colormap(hot);
colormap(cmap);
caxis([-40 40])
ax = gca;
ax.XTick = [];
ax.YTick = [];
ax.FontSize = 20;

pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig1/figs/selectivity/';
fn = 'selslice';
mysavefig(f,pth,fn);

%%

mask = obj.time > 0 & obj.time < 0.5;
means = mean(sel(mask,:),1);
[~,ix] = sort(means);
selsorted = sel(:,ix);

f = figure;
imagesc(obj.time,1:size(selsorted,2),selsorted')

% cmap = flipud(colormap(hot));
cmap = colormap(hot);

colormap(cmap)
c = colorbar;
caxis([-40 40])
ax = gca;
ax.FontSize = 20;
c.Label.String = 'Selectivity (spks/sec)';
xlabel('Time (s) from go cue')
ylabel('Neurons')










