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

psth = normalizePSTH(objs{1});
% psth = objs{1}.psth;
ev.sample = objs{1}.bp.ev.sample;
ev.delay = objs{1}.bp.ev.delay;
ev.goCue = objs{1}.bp.ev.goCue;
ev.(params.alignEvent) = objs{1}.bp.ev.(params.alignEvent);
for sessix = 2:numel(objs)
    obj = objs{sessix};
    
    temppsth = normalizePSTH(obj);
%     temppsth = obj.psth;
    
    psth = cat(2,psth,temppsth);
    
    ev.sample = [ev.sample ; obj.bp.ev.sample];
    ev.delay = [ev.delay ; obj.bp.ev.delay];
    ev.goCue = [ev.goCue ; obj.bp.ev.goCue];
    ev.(params.alignEvent) = [ev.(params.alignEvent) ; obj.bp.ev.(params.alignEvent)];
    
end

psth(isnan(psth)) = 0;
psth(isinf(psth)) = 0;

%% ACTIVITY MODES

rez.time = objs{1}.time;
rez.psth = psth;
rez.condition = params.condition;
rez.alignEvent = params.alignEvent;
rez.ev = ev;

%% cd early

% TODO: write calcCD_multi. Should return mu and sd for each loaded obj -
% concatenate them. 
% then, calculate each mode
% Then concatenate PSTHs, then project and plot

% cond{1} = params.modecondition{1};
% cond{2} = params.modecondition{2};
% epoch = 'early';
% 
% rez.cdEarly_mode = calcCD_multi(objs,params,cond,epoch);
% 
% clear cond

%% cd late mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
epoch = 'latedelay';

rez.cdLate_mode = calcCD_multi(objs,params,cond,epoch);

clear cond

%% cd go mode
cond{1} = params.modecondition{1};
cond{2} = params.modecondition{2};
epoch = 'go';

rez.cdGo_mode = calcCD_multi(objs,params,cond,epoch);

clear cond

%% orthogonalize

[fns,~] = patternMatchCellArray(fieldnames(rez),{'mode'},'all');
modes = zeros(numel(rez.cdLate_mode),numel(fns));
for i = 1:numel(fns)
    modes(:,i) = rez.(fns{i});
end

orthModes = gschmidt(modes);

for i = 1:numel(fns)
    rez.(fns{i}) = orthModes(:,i);
end


%% PLOTS

clrs = getColors();


% plot correct trials alone
plt.title = 'Correct Trials';
plt.legend = {'Right Hit','Left Hit'};
plt.conditions = [1,2];
plt.lw = [4 4];
plt.smooth = 51;
plt.colors = {clrs.rhit,clrs.lhit};
plt.save = 0;
plotAllCDs(rez, rez.ev, params.alignEvent, plt) 

% % plot correct trials and error trials
% plt.title = 'Correct and Error Trials';
% plt.legend = {'Right Hit','Left Hit','Right Error', 'Left Error'};
% plt.conditions = [1,2,5,6];
% plt.lw = [4 4 3 3];
% plt.smooth = 51;
% plt.colors = {clrs.rhit,clrs.lhit,clrs.rmiss,clrs.lmiss};
% plotAllCDs(rez, obj.bp.ev, params.alignEvent, plt) 
% 



%% selectivity

selpsth = objs{1}.psth;
for sessix = 2:numel(objs)
    obj = objs{sessix};    
    selpsth = cat(2,selpsth,obj.psth);
end

selpsth(isnan(psth)) = 0;
selpsth(isinf(psth)) = 0;


cond = [1,2];
latent_cdlate = zeros(size(selpsth,1),numel(cond));
latent_cdgo = zeros(size(selpsth,1),numel(cond));
for i = 1:numel(cond)
    c = cond(i);
    latent_cdlate(:,i) = selpsth(:,:,c)*rez.cdLate_mode;
    latent_cdgo(:,i) = selpsth(:,:,c)*rez.cdGo_mode;
end

sm = 31;

psth_selectivity = mySmooth(selpsth(:,:,1) - selpsth(:,:,2),sm);
cdlate_selectivity = mySmooth(latent_cdlate(:,1) - latent_cdlate(:,2),sm);
cdgo_selectivity = mySmooth(latent_cdgo(:,1) - latent_cdgo(:,2),sm);


sample = mode(rez.ev.sample) - mode(rez.ev.(params.alignEvent));
delay  = mode(rez.ev.delay) - mode(rez.ev.(params.alignEvent));

lw = 4;
figure; 
plot(rez.time,mean(psth_selectivity.^2,2),'k','LineWidth',lw)
hold on
plot(rez.time,mean(cdlate_selectivity.^2,2),'g','LineWidth',lw)
plot(rez.time,mean(cdgo_selectivity.^2,2),'m','LineWidth',lw)
plot(rez.time,(cdgo_selectivity.^2 + cdlate_selectivity.^2),'c','LineWidth',lw)

xline(sample,'k--','LineWidth',0.5);
xline(delay,'k--','LineWidth',0.5);
xline(0,'k--','LineWidth',0.5);

xlabel('Time (s) from go cue')
ylabel('Squared Selectivity')
legend('Total selectivity','late','go','late + go')
xlim([rez.time(1)+0.2,rez.time(end)])
ax = gca;
ax.FontSize = 20;


%% correlation pop selectivity vector

corr_matrix_selectivity = zeros(size(psth_selectivity,1),size(psth_selectivity,1));

for i = 1:size(corr_matrix_selectivity,1)
    for j = 1:size(corr_matrix_selectivity,1)
        temp = corrcoef(psth_selectivity(i,:),psth_selectivity(j,:));
        corr_matrix_selectivity(i,j) = temp(1,2);
    end
end

%% plot selectivity correlation matrix

figure; hold on;
imagesc(obj.time,obj.time,corr_matrix_selectivity);
colorbar; caxis([0 max(max(corr_matrix_selectivity))]);

lw = 4;
xline(sample,'w--','LineWidth',lw); yline(sample,'w--','LineWidth',lw)
xline(delay,'w--','LineWidth',lw); yline(delay,'w--','LineWidth',lw)
xline(0,'w--','LineWidth',lw); yline(0,'w--','LineWidth',lw)

xlim([rez.time(1)+0.2,rez.time(end)]);
ylim([rez.time(1)+0.2,rez.time(end)])

ax = gca;
ax.FontSize = 20;
hold off
colormap(hot)



