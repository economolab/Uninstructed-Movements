% Finding CDContext from neural activity that resides within the Null and Potent spaces
% Then finding selectivity between CDContext from full neural pop and
% CDContext found from null/potent reconstructions
clear,clc,close all

whichcomp = 'LabPC';                                                % LabPC or Laptop

% Base path for code depending on laptop or lab PC
if strcmp(whichcomp,'LabPC')
    basepth = 'C:\Users\Jackie Birnbaum\Documents\Code';
elseif strcmp(whichcomp,'Laptop')
    basepth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code';
end

% add paths
utilspth = [basepth '\Munib Uninstruct Move\uninstructedMovements_v2'];
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'figNP')));
figpth = [basepth  '\Uninstructed-Movements\Fig 3'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context switching')));
figpth = [basepth  '\Uninstructed-Movements\Fig 6'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context_funcs')));
figpth = [basepth  '\Uninstructed-Movements\Fig 5'];
addpath(genpath(fullfile(figpth,'funcs')));
figpth = [basepth  '\Uninstructed-Movements\Fig 2'];
addpath(genpath(fullfile(figpth,'funcs')));
%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'hit&~stim.enable&~autowater'};               % all 2AFC hits, no stim
params.condition(end+1) = {'hit&~stim.enable&autowater'};                % all AW hits, no stim
params.condition(end+1) = {'miss&~stim.enable&~autowater'};              % error 2AFC, no stim, aw off
params.condition(end+1) = {'miss&~stim.enable&autowater'};               % error AW, no stim

params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % all 2AFC hits, ~early, no stim
params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};                % all AW hits, ~early,no stim

params.tmin = -3;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;

params.bctype = 'reflect'; % options are : reflect  zeropad  none
%% SPECIFY DATA TO LOAD

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

meta = [];

% --- ALM --- 
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written
%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);
%%
% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end
%% Find Null and Potent Spaces using all trials (no train/test split yet)
clearvars -except obj meta params me sav

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -----------------------------------------------------------------------
for sessix = 1:numel(meta)
    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix));
    zscored(sessix).trialdat =  trialdat_zscored;

    % -- Calculate the null and potent spaces for each session
    cond2use = [2 3 4 5];   % All 2AFC hit/miss trials, all AW hit/miss trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
    nullalltime = 0;        % use all time points to estimate null space if 1
    AWonly = 0;             % use only AW to find null and potent spaces 
    delayOnly = 0;          % use only delay period to find null and potent spaces
    cond2proj = [6 7];      % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj,nullalltime,AWonly,delayOnly); 
end
%% Split trials into Move vs. Non-Move (during presample)

% Specify times in which you want to find move/non-move trials
times.trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
times.startix = find(obj(1).time>times.trialstart,1,'first');
times.samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
times.stopix = find(obj(1).time<times.samp,1,'last');

cond2use = [6 7];                                                           % Conditions (in reference to params.trialid) that you want to find move and non-move trials for
trials2cutoff = 40;                                                         % Trials to cut-off at the end of the session
% MoveNonMove = (1 x nSessions) struct with fields 'noMove' and 'Move'
% Each field has subfields 'afc' and 'aw'
% 'afc' and 'aw' are (n x 1) where 'n' = the number of move or non-move trials in that context
MoveNonMove = findMoveNonMoveTrix(meta, obj, trials2cutoff, cond2use, params, me, times);
  %% Get train/test split for data
clearvars -except obj meta params rez zscored me MoveNonMove times

% Store indices for train and test trials in variable called 'testsplit'
trainPct = 0.5;     % Percentage of trials being used for train vs test
condfns = {'afc','aw'};
movefns = {'noMove','Move','all'};

testsplit = getTestTrials_MnM(condfns,movefns,trainPct,MoveNonMove);

%testsplit = getTestTrials(params,cond2use,trainPct);
%% Find CDContext using training data (all trials, not separated by move/non-move)
% Projections are only for test data (all trials and also separated by Move/NonMove)
condfns = {'afc','aw'};
popfns = {'null','potent','fullpop'};
movefns = {'noMove','Move','all'};

% New field added to 'cd_null' = 'testsingleproj'
% testsingleproj has fields 'noMove', 'Move', 'all' indicating whether these single trial projections are from noMove trials, move trials, or all trials
% Each field then has a (1 x 2) cell array.  First cell contains single test trial projections from the 2AFC context. 
% Second cell contains projections from AW context
[cd_null, cd_potent, cd_context] = CDContext_AllSpaces_MnM(obj,meta,rez,popfns,condfns,movefns,testsplit,params);

clearvars -except cd_context cd_null cd_potent obj params meta me rez zscored testsplit nSplits times
%% Reorganize into a simpler structure
% Grouped = (1 x nSessions) struct with fields 'fullpop', 'null', 'potent'
% Each field has subfields 'noMove','Move','all'
% 'noMove','Move', and 'all' are (1 x 2) cell. Where first cell contains single trial projections onto CDContext for 2AFC trials (of that move condition)
% and second cell contains for AW trials

popfns = {'fullpop','null','potent'};
movefns = {'noMove','Move','all'};

grouped = reorganizeMnM(meta, popfns, movefns, cd_context, cd_null, cd_potent);

clearvars -except cd_context cd_null cd_potent obj params meta me rez zscored testsplit times grouped popfns movefns
%% Group across all sessions
% Find average CDContext in all move conditions and task conditions
% Find selectivity (btw 2AFC and AW)

% all_grouped = (1 x nSessions) struct with fields 'fullpop', 'null', 'potent'
% Each field has subfields 'noMove','Move','all'
% 'noMove','Move', and 'all' are (1 x 2) cell. Where first cell contains single trial projections onto CDContext for 2AFC trials (of that move condition)
% and second cell contains for AW trials
sm=50;
all_grouped = combineSessions_grouped(meta,grouped,sm, popfns, movefns);
%% Plot average CDContext across all sessions for each context

colors = getColors();
alph = 0.2;             % Shading opacity for error bars
%%% OLD VERSION OF FIGURE %%%
% LinePlot_CDGrouped_MoveNonMove(meta,ngroups,all_grouped,trialstart,samp,alph,colors,obj)
%% Plot average selectivity in CDContext across all sessions for each context
% figure();
% LinePlot_SelGrouped_MoveNonMove_V1(meta,ngroups,all_grouped,trialstart,samp,alph,colors,obj)

figure();
LinePlot_SelGrouped_MoveNonMove_V2(meta,all_grouped,times,alph,colors,obj,movefns,popfns)
%% Get the average presample selectivity in CDCont for move and non-move trials
for ii = 1:length(popfns)
    cont = popfns{ii};
    for gg = 1:length(movefns)
        presampavg.(cont).(movefns{gg}) = mean(all_grouped.(cont).(movefns{gg}).selectivity(times.startix:times.stopix,:),1,'omitnan');
    end
end
%% Normalize the selectivity to the max presample selectivity in each "condition" i.e. fullpop, null, potent
for ii = 1:length(popfns)
    cont = popfns{ii};
    temp = [presampavg.(cont).noMove, presampavg.(cont).Move];              % Take all presamp average values for CDCont, in move vs non-move
    maxsel = max(abs(temp));                                                % Take the maximum of the absolute value of these
   
    presampavgnorm.(cont).noMove = abs(presampavg.(cont).noMove)./maxsel;   % Normalize all presamp avgs from this condition to this value 
    presampavgnorm.(cont).Move = abs(presampavg.(cont).Move)./maxsel;       % For Move as well
end
%% Do all t-tests (paired)
sigcutoff = 0.01;
allhyp = [];            % (3 x 1).  First row = full pop.  Second row = null. Third row = potent.             
for ii = 1:length(popfns)
    cont = popfns{ii};
    hyp.(cont) = ttest(presampavgnorm.(cont).noMove,presampavgnorm.(cont).Move,'Alpha',sigcutoff);
end
%% Bar plot + scatter plot of avg presample CDContext on move trials vs non-move trials
figure();
plotBarPlot_Scatter_WLines(presampavgnorm)
%% Plotting functions

function LinePlot_CDGrouped_MoveNonMove(meta,ngroups,all_grouped,trialstart,samp,alph,colors,obj)
nSessions = length(meta);
cnt = 1;
for ii = 1:3
    if ii==1
        cont = 'fullpop';
        yl = [-5 12];
    elseif ii==2
        cont = 'null';
        yl = [-0.1 0.15];
    elseif ii==3
        cont = 'potent';
        yl = [-0.4 0.1];
    end
    for gg = 1:ngroups
        if gg==1
            movement = 'Move trials';
        else
            movement = 'Non-Move trials';
        end
        for cc = 1:2
            if cc == 1
                trialcont = 'afc';
                col = colors.afc;
            else
                trialcont = 'aw';
                col = colors.aw;
            end
            subplot(3,ngroups,cnt)
            ax = gca;
            toplot = mean(all_grouped.(cont).(trialcont){gg},2,'omitnan');
            err = std(all_grouped.(cont).(trialcont){gg},0,2,'omitnan')./sqrt(nSessions);
            %err = 1.96*(std(all_grouped.(cont).(trialcont){gg},0,2,'omitnan')./sqrt(nSessions));
            shadedErrorBar(obj(1).time,toplot,err,{'Color',col,'LineWidth',2},alph,ax); hold on;
        end
        title([cont ';  ' movement])
        %ylim(yl)
        if ii~=1
            set(ax, 'YDir','reverse')
        end
        xlim([trialstart 2.5])
        xline(0,'k--','LineWidth',1)
        xline(samp,'k--','LineWidth',1)
        cnt = cnt+1;
    end
end
end

function LinePlot_SelGrouped_MoveNonMove_V1(meta,ngroups,all_grouped,trialstart,samp,alph,colors,obj)
nSessions = length(meta);
cnt = 1;
for ii = 1:3
    if ii==1
        cont = 'fullpop';
        yl = [-5 12];
        col = [0.25 0.25 0.25];
    elseif ii==2
        cont = 'null';
        yl = [-0.15 0.15];
        col = colors.null;
    elseif ii==3
        cont = 'potent';
        yl = [-0.45 0.1];
        col = colors.potent;
    end
    for gg = 1:ngroups
        if gg==1
            movement = 'Move trials';
        else
            movement = 'Non-Move trials';
        end

        subplot(3,ngroups,cnt)
        ax = gca;
        toplot = mean(all_grouped.(cont).selectivity{gg},2,'omitnan');
        err = std(all_grouped.(cont).selectivity{gg},0,2,'omitnan')./sqrt(nSessions);
        %err = 1.96*(std(all_grouped.(cont).(trialcont){gg},0,2,'omitnan')./sqrt(nSessions));
        shadedErrorBar(obj(1).time,toplot,err,{'Color',col,'LineWidth',2},alph,ax);

        title([cont '; ' movement])
%         ylim(yl)
        if ii~=1
            set(ax, 'YDir','reverse')
        end
        xlim([trialstart 2.5])
        xline(0,'k--','LineWidth',1)
        xline(samp,'k--','LineWidth',1)
        cnt = cnt+1;
    end
end
end

function LinePlot_SelGrouped_MoveNonMove_V2(meta,all_grouped,times,alph,colors,obj,movefns,popfns)
nSessions = length(meta);
cnt = 1;
for po = 1:length(popfns)
    cont = popfns{po};
    switch cont
        case 'fullpop'
        yl = [0 8.5];
        col = [0.25 0.25 0.25];
        case 'null'
        yl = [0 0.35];
        col = colors.null;
        case 'potent'
        yl = [0 0.35];
        col = colors.potent;
    end
    subplot(3,2,cnt)
    ax = gca;
    toplot = mean(all_grouped.(cont).all.selectivity,2,'omitnan');
    err = std(all_grouped.(cont).all.selectivity,0,2,'omitnan')./sqrt(nSessions);
    %err = std(all_grouped.(cont).all.selectivity,0,2,'omitnan')./sqrt(nSessions);
    shadedErrorBar(obj(1).time,toplot,err,{'Color',col,'LineWidth',2},alph,ax);

    ylim(yl)
    xline(times.samp,'k--','LineWidth',1)
    xlim([times.trialstart 0])
    ylabel(cont)
    title('All trials')
    cnt = cnt+1;

    for gg = 1:2
        if gg==2                            % Move trials
            style = '-';
            alp = alph;
        elseif gg==1                        % Non-Move trials
            style = '--';
            alp = alph-0.1;
        end

        subplot(3,2,cnt)
        ax = gca;
        toplot = mean(all_grouped.(cont).(movefns{gg}).selectivity,2,'omitnan');
        err = std(all_grouped.(cont).(movefns{gg}).selectivity,0,2,'omitnan')./sqrt(nSessions);
        %err = 1.96*(std(all_grouped.(cont).selectivity{gg},0,2,'omitnan')./sqrt(nSessions));
        shadedErrorBar(obj(1).time,toplot,err,{'Color',col,'LineWidth',2,'LineStyle',style},alp,ax);
        hold on;

        ylim(yl)
%         if ii~=1
%             set(ax, 'YDir','reverse')
%         end
%        xline(0,'k--','LineWidth',1)
        xline(times.samp,'k--','LineWidth',1)
        xlim([times.trialstart 0])
        ylabel('Selectivity (a.u.)')
        legend({movefns{1} movefns{2}})
    end
    cnt = cnt+1;
end
end

function plotBarPlot_Scatter_WLines(presampavgnorm)
X = [1,2,4,5,7,8];
row1 = []; row2 = []; row3 =[];
movefns = {'Move','noMove'};
for gg = 1:length(movefns)
    row1 = [row1,mean(presampavgnorm.fullpop.(movefns{gg}))];
    row2 = [row2,mean(presampavgnorm.null.(movefns{gg}))];
    row3 = [row3,mean(presampavgnorm.potent.(movefns{gg}))];
end
y = [row1, row2, row3];
figure();
bar(X,y);
hold on;
cnt = 1;
for ii  = 1:3
    switch ii
        case 1
            cont = 'fullpop';
            xx = [1,2];
        case 2
            cont = 'null';
            xx = [4,5];
        case 3
            cont = 'potent';
            xx= [7,8];
    end
    for gg = 1:length(movefns)
        scatter(X(cnt),presampavgnorm.(cont).(movefns{gg}),25,[0 0 0],'filled','MarkerEdgeColor','black')
        cnt = cnt+1;
    end
    for sessix = 1:length(meta)
        plot(xx,[presampavgnorm.(cont).Move(sessix),presampavgnorm.(cont).noMove(sessix)],'Color','black')
    end
end
end
%%
function grouped = reorganizeMnM(meta, popfns, movefns, cd_context, cd_null, cd_potent)
for sessix = 1:length(meta)
    for po = 1:length(popfns)
        fn = popfns{po};
        switch fn
            case 'fullpop'
                cd = cd_context(sessix).testsingleproj;
            case 'null'
                cd = cd_null(sessix).testsingleproj;
            case 'potent'
                cd = cd_potent(sessix).testsingleproj;
        end
        for mo = 1:length(movefns)
            grouped(sessix).(fn).(movefns{mo}) = cd.(movefns{mo});
        end
    end
end
end


function all_grouped = combineSessions_grouped(meta,grouped,sm, popfns, movefns)
for ii = 1:length(popfns)
    cont = popfns{ii};
    for mo = 1:length(movefns)
        temp = [];
        for sessix = 1:length(meta)
            tempafc = mySmooth(mean(grouped(sessix).(cont).(movefns{mo}){1},2,'omitnan'),sm);
            tempaw = mySmooth(mean(grouped(sessix).(cont).(movefns{mo}){2},2,'omitnan'),sm);
            sel = tempafc-tempaw;
            temp = [temp,sel];
        end
        all_grouped.(cont).(movefns{mo}).selectivity = temp;
            
        for cc = 1:2
            if cc == 1
                trialcont = 'afc';
            else
                trialcont = 'aw';
            end
            temp = [];
            for sessix = 1:length(meta)
                temp = [temp,mean(grouped(sessix).(cont).(movefns{mo}){cc},2,'omitnan')];
            end
            all_grouped.(cont).(movefns{mo}).(trialcont) = temp;
        end
    end
end
end

