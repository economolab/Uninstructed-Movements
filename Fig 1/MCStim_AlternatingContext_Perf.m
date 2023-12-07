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
params.condition(1)     = {'(hit|miss|no)'};                             % all trials

params.condition(end+1) = {'L&~stim.enable&~autowater&~early'};             % right, no stim, 2AFC
params.condition(end+1) = {'L&stim.enable&~autowater&~early'};              % right, stim, 2AFC
params.condition(end+1) = {'R&~stim.enable&~autowater&~early'};             % left, no stim, 2AFC
params.condition(end+1) = {'R&stim.enable&~autowater&~early'};              % left, stim, 2AFC

params.condition(end+1) = {'L&~stim.enable&autowater&~early'};              % right, no stim, AW
params.condition(end+1) = {'L&stim.enable&autowater&~early'};               % right, stim, AW
params.condition(end+1) = {'R&~stim.enable&autowater&~early'};              % left, no stim, AW
params.condition(end+1) = {'R&stim.enable&autowater&~early'};               % left, stim, AW
%% SPECIFY DATA TO LOAD

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end
meta = [];
meta = loadMAH13_MCStim(meta,datapth);
meta = loadMAH14_MCStim(meta,datapth);
meta = loadMAH20_MCStim(meta,datapth);
meta = loadMAH21_MCStim(meta,datapth);
%% LOAD DATA

% ----------------------------------------------
% -- Behavioral and Video Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
disp('Loading behavior objects')
[obj,params] = loadBehavSessionData(meta,params);
%% Find trials with a correct lick within a certain latency from the go cue
conds2use = 2:9;                          % L,R control DR; L,R stim DR; L,R control WC; L,R stim WC
lickCutoff = 0.6;                         % Window after go cue that you are looking for correct licks
trialCutoff = 40;                         % How many trials you want to cut off of the end of the session (to account for loss of motivation/engagement)
blah = [];
for sessix = 1:length(meta)               % For each session...
    lastTrial = obj(sessix).bp.Ntrials-trialCutoff;        % Get rid of the last 'trialCutoff' trials in the session
    for c = 1:length(conds2use)                            % For each condition...                                         
        cond = conds2use(c);
        condtrix = params(sessix).trialid{cond};              % Get the trials for that condition
        condtrix = condtrix(condtrix<lastTrial);              % Exclude trials that are beyond the trialCutoff
        lickInit = false(length(condtrix),1);                 % [# of trials in condition x 1]
        for t = 1:length(condtrix)                            % For every trial in condition...                             
            currtrial = condtrix(t);
            allLicks = [obj(sessix).bp.ev.lickL{currtrial},obj(sessix).bp.ev.lickR{currtrial}];     % Get the times of all licks in this trial
            allLicks = sort(allLicks,'ascend');                 % Sort the licks in chronological order
            if ~isempty(allLicks)                               % As long as there were licks in this trial...
                go = obj(sessix).bp.ev.goCue(currtrial);          % Find the time of the go cue on this trial
                postGoLicks = allLicks((allLicks>go));            % Find the licks that occurred after the go cue
                if ~isempty(postGoLicks)                              % As long as there were licks after the go cue on this trial...
                    firstLick_aligned = postGoLicks(1)-go;               % Get the time of the first lick that occurs after the go cue (time is aligned to the go cue)
                    if firstLick_aligned<lickCutoff && obj(sessix).bp.hit(currtrial)      % If that first lick occurred before the cutoff time ('lickCutoff') and it was a correct trial...
                        lickInit(t) = 1;                                                    % 0 = no correct lick before LickCutoff time; 1 = there is a correct lick before LickCutoff
                    end
                end
            end
        end
        perf(c) = sum(lickInit)/length(condtrix);              % Number of trials with a correct lick within the time window / all trials in that condition
    end
    rez(sessix).trixWLicks = perf;
end
%% Concatenate across sessions
clearvars -except obj meta rez params lickCutoff 
perf_all = [];
for sessix = 1:length(meta)
    perf_all = [perf_all;rez(sessix).trixWLicks];
end
%% INCLUSION CRITERIA: Omit sessions where 2AFC stim didn't work or where performance was bad
sess2omit = false(length(meta),1);
Rperf = NaN(length(meta),1);  Reff = NaN(length(meta),1);
Lperf = NaN(length(meta),1);  Leff = NaN(length(meta),1);

LAFCctrlCond = 1;
RAFCctrlCond = 3;
EffectCutoff = 0;
perfCutoff = 0.55;
for sessix = 1:length(meta)
    temp = perf_all(sessix,:);
    Reff(sessix) = temp(RAFCctrlCond)-temp(RAFCctrlCond+1);
    Leff(sessix) = temp(LAFCctrlCond)-temp(LAFCctrlCond+1);
    if Reff(sessix) < EffectCutoff ||  Leff(sessix)< EffectCutoff
        sess2omit(sessix) = 1;
    end
    bp = obj(sessix).bp;
    Rperf(sessix) = sum(bp.R&bp.hit&~bp.autowater)/sum(bp.R&~bp.autowater); 
    Lperf(sessix) = sum(bp.L&bp.hit&~bp.autowater)/sum(bp.L&~bp.autowater);
    if Rperf(sessix)<perfCutoff || Lperf(sessix)<perfCutoff
        sess2omit(sessix) = 1;
    end
end
perf_all(sess2omit,:) = [];

% Do t-test
clearvars -except obj meta rez params perf_all lickCutoff sess2omit
%perf_all = (sessions x conditions)
perf_all = 100*perf_all;
%% Plot
anmNames = {'MAH13','MAH13','MAH13','MAH13','MAH13','MAH13',...
    'MAH14','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14','MAH14',...
    'MAH20','MAH20','MAH20','MAH20','MAH20','MAH20','MAH20',...
    'MAH21','MAH21','MAH21','MAH21','MAH21','MAH21','MAH21','MAH21',};
anmNames(sess2omit) = [];
cols = getColors();
sigcutoff = 0.05;
[AFCpval,AWpval] = plotPerfV1(cols,perf_all,lickCutoff,anmNames,sigcutoff);

%plotPerfV2(colors,perf_all,lickCutoff,anmNames)
%% Print summary statistics 
disp('---Summary statistics for MC go cue photoinhibition---')
disp(['p-values for 2AFC t-tests -- L: ' num2str(AFCpval(1)) ' ; R: ' num2str(AFCpval(2))])
disp(['p-values for AW t-tests -- L: ' num2str(AWpval(1)) ' ; R: ' num2str(AWpval(2))])
disp('Paired t-test')
t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
disp(t)
%%
function [AFCpval,AWpval] = plotPerfV1(cols,perf_all,lickCutoff,anmNames,sigcutoff)
figure();
subplot(1,2,1)
uniqueAnm = unique(anmNames);
nAnimals = length(uniqueAnm);

conds2plot = 1:4;
for x = conds2plot
    switch x
        case 1
            facecol = cols.lhit;
            edgecol = [1 1 1];
        case 2
            facecol = [1 1 1];
            edgecol = cols.lhit;
        case 3
            facecol = cols.rhit;
            edgecol = [1 1 1];
        case 4
            facecol = [1 1 1];
            edgecol = cols.rhit;
    end

    bar(x,mean(perf_all(:,x)),'FaceColor',facecol,'EdgeColor',edgecol); hold on;
    for anm = 1:nAnimals
        curranm = uniqueAnm{anm};
        anmix = find(strcmp(anmNames,curranm));
        switch anm
            case 1
                shape = 'o';
            case 2
                shape = '^';
        end
        xx = x*ones(length(anmix),1);
        scatter(xx,perf_all(anmix,x),'filled',shape,'MarkerFaceColor','black');
    end
end
for sessix = 1:size(perf_all,1)
    plot(1:2,perf_all(sessix,1:2),'Color','black')
    plot(3:4,perf_all(sessix,3:4),'Color','black')
end
for test = 1:(length(conds2plot)/2)
    x = perf_all(:,(test*2)-1);
    y = perf_all(:,test*2);
    [hyp,AFCpval(test)] = ttest(x,y,'Alpha',sigcutoff);
    disp(num2str(hyp))
    if hyp&&test==1
        scatter(1.5,105,30,'*','MarkerEdgeColor','black')
    elseif hyp&&test==2
        scatter(3.5,105,30,'*','MarkerEdgeColor','black')
    end
end
xlim([0 5])
xticks([1,2,3,4])
xticklabels({'L ctrl', 'L stim','R ctrl', 'R stim'})
ylabel(['Proportion of trials w/ lick within ' num2str(lickCutoff) ' (s) of goCue'])
title('2AFC')

subplot(1,2,2)
conds2plot = 5:8;
for x = 1:length(conds2plot)
    cond = conds2plot(x);

    switch x
        case 1
            facecol = cols.lhit_aw;
            edgecol = [1 1 1];
        case 2
            facecol = [1 1 1];
            edgecol = cols.lhit_aw;
        case 3
            facecol = cols.rhit_aw;
            edgecol = [1 1 1];
        case 4
            facecol = [1 1 1];
            edgecol = cols.rhit_aw;
    end

    bar(x,mean(perf_all(:,cond)),'FaceColor',facecol,'EdgeColor',edgecol); hold on;
    for anm = 1:nAnimals
        curranm = uniqueAnm{anm};
        anmix = find(strcmp(anmNames,curranm));
        switch anm
            case 1
                shape = 'o';
            case 2
                shape = '^';
        end
        xx = x*ones(length(anmix),1);
        scatter(xx,perf_all(anmix,cond),'filled',shape,'MarkerFaceColor','black');
    end
end
for sessix = 1:size(perf_all,1)
    plot(1:2,perf_all(sessix,conds2plot(1:2)),'Color','black')
    plot(3:4,perf_all(sessix,conds2plot(3:4)),'Color','black')
end
for test = 1:(length(conds2plot)/2)
    x = perf_all(:,conds2plot((test*2)-1));
    y = perf_all(:,conds2plot(test*2));
    [hyp,AWpval(test)] = ttest(x,y,'Alpha',sigcutoff);
    disp(num2str(hyp))
    if hyp&&test==1
        scatter(1.5,105,30,'*','MarkerEdgeColor','black')
    elseif hyp&&test==2
        scatter(3.5,105,30,'*','MarkerEdgeColor','black')
    end
end
xlim([0 5])
xticks([1,2,3,4])
xticklabels({'L ctrl', 'L stim','R ctrl', 'R stim'})
ylabel(['Proportion of trials w/ lick within ' num2str(lickCutoff) ' (s) of waterDrop'])
title('Autowater')
end


