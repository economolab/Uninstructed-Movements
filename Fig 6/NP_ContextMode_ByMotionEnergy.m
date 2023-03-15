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
addpath(genpath(fullfile(utilspth,'fig3')));
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
%% Null and Potent Space
clearvars -except obj meta params me sav

% -----------------------------------------------------------------------
% -- Curate Input Data --
% zscore single trial neural data (time*trials,neurons), for all trials
% -- Calculate null and potent spaces --
% null space from quiet time points
% potent space from moving time points
% -- Calculate coding directions from null and potent spaces --
% -- Calculate coding directions from full neural population --
% -----------------------------------------------------------------------

for sessix = 1:numel(meta)
    % -- input data
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat, obj(sessix));

    % -- Calculate the null and potent spaces for each session
    cond2use = [2 3 4 5];   % All 2AFC hit/miss trials, all AW hit/miss trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
    nullalltime = 0;        % use all time points to estimate null space if 1
    AWonly = 0;             % use only AW to find null and potent spaces 
    delayOnly = 0;          % use only delay period to find null and potent spaces
    cond2proj = [6 7];      % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj,nullalltime,AWonly,delayOnly);

    % -- Find coding dimensions from RECONSTRUCTED full neural activity which is reconstructed from the null and potent spaces
    cond2use = [1 2];            % (NUMBERING ACCORDING TO THE CONDITIONS PROJECTED INTO NULL AND POTENT SPACES, i.e. which of the conditions specified in 'cond2proj' above do you want to use?)
    cond2proj = [1 2];           % 2AFC hits/misses, AW hits/misses(corresponding to null/potent psths in rez)
    cond2use_trialdat = [6 7];   % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    cd_null(sessix) = getCodingDimensions_Context(rez(sessix).recon_psth.null,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);
    cd_potent(sessix) = getCodingDimensions_Context(rez(sessix).recon_psth.potent,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);

    % Calc CDContext from full neural pop
    cd_context(sessix) = getCodingDimensions_Context(obj(sessix).psth,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);
end
%% Get single trial projections onto CDContext
orthogonalize = 'non-orthog';                                       % Set to orthogonalize if you want the projections to be onto the orthogonalized CDs
disp('----Projecting single trials onto CDContext----')
cd = 'context';
cd_context = getSingleTrialProjs(cd_context,obj,cd,orthogonalize);
%% Get single trial NP projections onto CDContext
cd = 'context';
[cd_null,cd_potent] = getNPSingleTrialProjs(obj,cd,cd_null,cd_potent,rez);
%% Get trials with a lot of presample motion energy and little presample motion energy
clearvars -except obj meta params me sav rez cd_null cd_potent cd_context

% Find the times corresponding to trial start and the sample period
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
startix = find(obj(1).time>trialstart,1,'first');
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
stopix = find(obj(1).time<samp,1,'last');

trials2cutoff = 40;                 % Trials to discount at the end of the session (for motion energy)
cond2use = [6,7];                   % 2AFC trials and AW trials
ngroups = 3;

grouped = groupCDbyME(meta, obj, trials2cutoff, cond2use, startix, stopix,cd_context,cd_null,cd_potent,ngroups,me);
%% Plot for each session
colors = getColors_Updated();
% for sessix = 1:length(meta)
%     cnt = 1;
%     for ii = 1:3
%         if ii==1
%             cont = 'fullpop';
%         elseif ii==2
%             cont = 'null';
%         elseif ii==3
%             cont = 'potent';
%         end
%         for gg = 1:ngroups
%             for cc = 1:2
%                 if cc == 1
%                     trialcont = 'afc';
%                     col = colors.afc;
%                 else
%                     trialcont = 'aw';
%                     col = colors.aw;
%                 end
%                 subplot(3,ngroups,cnt)
%                 plot(obj(1).time,grouped(sessix).(cont).(trialcont){gg},'Color',col); hold on;
%             end
%             legend('2AFC','AW')
%             title([cont '; Group ' num2str(gg)])
%             cnt = cnt+1;
%         end
%     end
%     pause
%     close all;
% end
%% Group across all sessions
sm=31;
ngroups = 2;            % Move vs non-move
all_grouped = combineSessions_grouped(ngroups,meta,grouped,sm);
%% Plot averages across all sessions
alph = 0.2;
LinePlot_CDGrouped_MoveNonMove(meta,ngroups,all_grouped,trialstart,samp,alph,colors,obj)
%%
LinePlot_SelGrouped_MoveNonMove(meta,ngroups,all_grouped,trialstart,samp,alph,colors,obj)
%% Get the average presample selectivity in CDCont for move and non-move trials
for ii = 1:3
    if ii==1
        cont = 'fullpop';
    elseif ii==2
        cont = 'null';
    elseif ii==3
        cont = 'potent';
    end
    for gg = 1:ngroups
        presampavg.(cont){gg} = mean(all_grouped.(cont).selectivity{gg}(startix:stopix,:),1,'omitnan');
    end
end
%% Normalize the selectivity to the max presample selectivity in each "condition" i.e. fullpop, null, potent
for ii = 1:3
    if ii==1
        cont = 'fullpop';
    elseif ii==2
        cont = 'null';
    elseif ii==3
        cont = 'potent';
    end
    temp = [presampavg.(cont){1}, presampavg.(cont){2}];
    maxsel = max(abs(temp));
    for gg = 1:ngroups
        presampavgnorm.(cont){gg} = abs(presampavg.(cont){gg})./maxsel;
    end
end

%% Do all t-tests (paired)
sigcutoff = 0.05;
allhyp = [];            % (3 x 1).  First row = full pop.  Second row = null. Third row = potent.  
                        
for ii = 1:3
    switch ii
        case 1
            cont = 'fullpop';
        case 2
            cont = 'null';
        case 3
            cont = 'potent';
    end
    hyp.(cont) = ttest(presampavgnorm.(cont){1},presampavgnorm.(cont){2},'Alpha',sigcutoff);
end
%%
X = [1,2,4,5,7,8];
row1 = []; row2 = []; row3 =[];
for gg = 1:ngroups
    row1 = [row1,mean(presampavgnorm.fullpop{gg})];
    row2 = [row2,mean(presampavgnorm.null{gg})];
    row3 = [row3,mean(presampavgnorm.potent{gg})];
end
y = [row1, row2, row3];
bar(X,y);
hold on;
cnt = 1;
for ii  = 1:3
    switch ii
        case 1
            cont = 'fullpop';
        case 2
            cont = 'null';
        case 3
            cont = 'potent';
    end
    for gg = 1:ngroups
        scatter(X(cnt),presampavgnorm.(cont){gg},25,[0 0 0],'filled','MarkerEdgeColor','black')
        cnt = cnt+1;
    end

end
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
            movement = 'Non-Move trials';
        else
            movement = 'Move trials';
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
        ylim(yl)
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

function LinePlot_SelGrouped_MoveNonMove(meta,ngroups,all_grouped,trialstart,samp,alph,colors,obj)
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
            movement = 'Non-Move trials';
        else
            movement = 'Move trials';
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
%%
function all_grouped = combineSessions_grouped(ngroups,meta,grouped,sm)
for ii = 1:3
    if ii==1
        cont = 'fullpop';
    elseif ii==2
        cont = 'null';
    elseif ii==3
        cont = 'potent';
    end
    for gg = 1:ngroups
        temp = [];
        for sessix = 1:length(meta)
            tempafc = mySmooth(grouped(sessix).(cont).afc{gg},sm);
            tempaw = mySmooth(grouped(sessix).(cont).aw{gg},sm);
            sel = tempafc-tempaw;
            temp = [temp,sel];
        end
        all_grouped.(cont).selectivity{gg} = temp;
            
        for cc = 1:2
            if cc == 1
                trialcont = 'afc';
            else
                trialcont = 'aw';
            end
            temp = [];
            for sessix = 1:length(meta)
                temp = [temp,grouped(sessix).(cont).(trialcont){gg}];
            end
            all_grouped.(cont).(trialcont){gg} = temp;
        end
    end
end
end

