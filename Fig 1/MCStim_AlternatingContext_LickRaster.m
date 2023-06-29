% Quantifying behavioral performance and cortical dependence in the
% Alternating Context Task
% -------------------------------------------------------------------------------------
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
figpth = [basepth  '\Uninstructed-Movements\Fig 1'];
addpath(genpath(fullfile(figpth,'Utils')));
%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% set conditions to calculate behavioral performance for 
params.condition(1) = {'L&~stim.enable&~autowater'};              % right, no stim, 2AFC
params.condition(end+1) = {'R&~stim.enable&~autowater'};              % left, no stim, 2AFC
params.condition(end+1) = {'L&stim.enable&~autowater'};               % right, stim, afc
params.condition(end+1) = {'R&stim.enable&~autowater'};               % left, stim, afc

params.condition(end+1) = {'L&~stim.enable&autowater'};              % right, no stim, aw
params.condition(end+1) = {'R&~stim.enable&autowater'};          % left, no stim, aw
params.condition(end+1) = {'L&stim.enable&autowater'};               % right, stim, aw
params.condition(end+1) = {'R&stim.enable&autowater'};               % left, stim, aw

params.alignEvent = 'goCue';
params.tmin = -2.5; 
params.tmax = 3;
params.dt = 1/200;

params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw'}};
params.smooth = 15;
params.advance_movement = 0;
params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this mu
%% SPECIFY DATA TO LOAD

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

meta = [];
meta = loadMAH13_MCStim(meta,datapth);
meta = loadMAH14_MCStim(meta,datapth);
%%
% ----------------------------------------------
% -- Behavioral and Video Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadBehavSessionData(meta,params);
for sessix = 1:length(obj)
    obj(sessix).time = params(sessix).tmin:params(sessix).dt:params(sessix).tmax;
end
%% Get kinematics
for sessix = 1:length(meta)
    kin(sessix) = getKinematics_NoME(obj(sessix),params(sessix));
end
%% 2AFC trials: top subplot = left trials (control on top, stim on bottom); bottom subplot = right trials (control on top, stim on bottom)
close all

sessix = 2; 

cols = getColors();
clrs.rhit = cols.rhit;
clrs.lhit = cols.lhit;

conds = [1 2 3 4];

figure(1)
plotLickRaster(sessix,clrs,obj,params,'2AFC',conds);
%% AW trials: top subplot = left trials (control on top, stim on bottom); bottom subplot = right trials (control on top, stim on bottom)
clrs.rhit = cols.rhit_aw;
clrs.lhit = cols.lhit_aw;

conds = [5 6 7 8];

figure(2)
plotLickRaster(sessix,clrs,obj,params,'AW',conds);
%% For supplement: show what tongue looks like during the stim period
sessix = 2;

trix2plot = 10;
smooth = 31;
offset = 5;

stim.stimstart = 0;
stim.stimstop = 1;
stim.stimepoch = 'Go cue';

kinfeat = 'tongue_length';
featix = find(strcmp(kin(sessix).featLeg,kinfeat));

cond2plot = [1 2 3 4];
condition = '2AFC';
figure();
plotKinTracking_CtrlvsStim(params,obj,sessix,cols,kin,offset,smooth,kinfeat,featix,trix2plot,cond2plot,stim,condition);

cond2plot = [5 6 7 8];
condition = 'AW';
figure();
plotKinTracking_CtrlvsStim(params,obj,sessix,cols,kin,offset,smooth,kinfeat,featix,trix2plot,cond2plot,stim,condition);
%% Calculate percentage of time where tongue is visible during stim period
cond2plot = 1:8;

stimepoch = 0;
stim.stimstart = 0;
stim.stimstop = 0.6;

for sessix = 1:length(meta)
    pctTime = NaN(1,length(cond2plot));

    for cond = 1:length(cond2plot)
        stimstart = stim.stimstart;
        startix = find(obj(sessix).time>stimstart,1,'first');
        stimstop = stim.stimstop;
        stopix = find(obj(sessix).time<stimstop,1,'last');

        condtrix = params(sessix).trialid{cond2plot(cond)};
        tempTime = NaN(1,length(condtrix));
        for t = 1:length(condtrix)
            currtrial = condtrix(t);
            featix = find(strcmp(kin(sessix).featLeg,'tongue_length'));
            tonguedat = kin(sessix).dat_std(startix:stopix,currtrial,featix);
            dur = length(tonguedat);
            tonguevis = length(find(~isnan(tonguedat)));
            tempTime(t) = tonguevis/dur;

            %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %         offset = 5;
            %         figure(cond)
            %         toplot = offset*t+tonguedat;
            %         if cond==1||cond==2
            %             plottime = 0.5+obj(sessix).time(startix:stopix);
            %         else
            %             plottime = 1+obj(sessix).time(startix:stopix);
            %         end
            %         plot(plottime,toplot); hold on
            %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        pctTime(cond) = mean(tempTime,'omitnan');
        %title(num2str(pctTime(cond)))
    end
    stimEffects(sessix).pctTime = pctTime;
end
%% Make bar plot comparing pct of time with tongue visible for each condition
pctTime = [];
for sess = 1:length(stimEffects)
    pctTime = [pctTime;stimEffects(sess).pctTime];
end

cols = getColors();
anmNames = {'MAH13','MAH13','MAH13','MAH13','MAH13','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14'};
sigcutoff = 0.05;

subplot(1,2,1)
condition = '2AFC';
starheight = 0.4;
temp = NaN(length(stimEffects),4);
temp(:,1) = pctTime(:,1);
temp(:,2) = pctTime(:,3);
temp(:,3) = pctTime(:,2);
temp(:,4) = pctTime(:,4);
AFCpvals = CtrlvsStimBarPlot(cols,temp,anmNames,sigcutoff,starheight,condition);
ylabel('Fraction of time')
ylim([0 0.42])
title('2AFC trials')
clear temp

subplot(1,2,2)
condition = 'AW';
starheight = 0.4;
temp = NaN(length(stimEffects),4);
temp(:,1) = pctTime(:,5);
temp(:,2) = pctTime(:,7);
temp(:,3) = pctTime(:,6);
temp(:,4) = pctTime(:,8);
AWpvals = CtrlvsStimBarPlot(cols,temp,anmNames,sigcutoff,starheight,condition);
ylabel('Fraction of time')
ylim([0 0.42])
title('Autowater trials')
clear temp
sgtitle(['Frac of time with tongue visible--from ' num2str(stim.stimstart) ' to ' num2str(stim.stimstop) ' (s)'])
%% Print summary statistics 
disp('---Pct of time with tongue visible for MC go cue photoinhibition---')
disp(['p-values for 2AFC t-tests -- L: ' num2str(AFCpvals(1)) ' ; R: ' num2str(AFCpvals(2))])
disp(['p-values for AW t-tests -- L: ' num2str(AWpvals(1)) ' ; R: ' num2str(AWpvals(2))])
disp('Paired t-test')
t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
disp(t)