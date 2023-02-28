% Find context selectivity across only Null and potent dimensions which are selective for context
% Take the difference between contexts for each selective dimension.  Take the mean
% across all selective dimensions within a session
% -------------------------------------------------------------------------------------
% Using all 2AFC and all AW trials to find the Null and Potent Spaces
% -------------------------------------------------------------------------------------
clear,clc,close all

% add paths
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig3')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 4';
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context switching')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 7';
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context_funcs')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 6';
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

params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};                   % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};                   % error left, no stim, aw off
params.condition(end+1) = {'~early&~stim.enable&~autowater'};                          %  no stim, aw off

params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};             % right hits, no stim, aw on
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};             % left hits, no stim, aw on
params.condition(end+1) = {'R&miss&~stim.enable&autowater'};                   % error right, no stim, aw on
params.condition(end+1) = {'L&miss&~stim.enable&autowater'};                   % error left, no stim, aw on
params.condition(end+1) = {'~early&~stim.enable&autowater'};                          %  no stim, aw on

params.condition(end+1) = {'hit&~early&~stim.enable&~autowater'};                          %  hit, no stim, aw off
params.condition(end+1) = {'hit&~early&~stim.enable&autowater'};                          %  hit, no stim, aw on

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

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

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
% -----------------------------------------------------------------------

for sessix = 1:numel(meta)
    % -- input data (the FR of each individual neuron will be converted to
    % a value between 0 and 1)
    trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix).trialdat, obj(sessix));

    % -- Calculate the null and potent spaces for each session
    cond2use = [2 3 7 8];   % All 2AFC hit/miss trials, all AW hit/miss trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
    nullalltime = 0;        % use all time points to estimate null space if 1
    AWonly = 0;             % use only AW to find null and potent spaces 
    delayOnly = 0;          % use only delay period to find null and potent spaces
    cond2proj = 2:11;       % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj,nullalltime,AWonly,delayOnly);
end
%%
% Find the times corresponding to trial start and the sample period
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>trialstart,1,'first');
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');
%% Find context selectivity in full neural population and in reconstructed neural pop

% Find cells that are significantly modulated by context in the presample period (ranksum test, p-value = 0.01)
% As in Inagaki et al., Cell, 2022 ('A midbrain-thalamus...')
% modulatedCells = (1 x nCells) array where 1 means the cell is context-selective and 0 means it is not
conds2comp = [12 13];

modulatedCells = [];
modparams.quals2excl = {'Poor','Multi','Noisy'};
modparams.sm = 30;                                                    % Amount that you want to smooth PSTHs by
modparams.measure = 'FR';                                             % Whether you want to compare firing rate ('FR') or spike counts ('spkCnt') 
modparams.subTrials = 35;
for sessix = 1:length(meta)                                 % For every session...   
    [selectivecells,spkdif] = findSelectiveCells(sessix,obj,modparams,meta,conds2comp);
    obj(sessix).selectiveCells = selectivecells;
    obj(sessix).spkdif = spkdif;
end
%% Calculate selectivity (spike rate difference) for all selective cells (on a given delay duration)
del2use = 0.9;
cond2use = [12,13];
smooth = 200;
for sessix = 1:length(meta)
    selectivity = findSelectivity_PSTH(sessix,obj,del2use,cond2use,smooth);
    obj(sessix).neuralselect = selectivity;
end
%% Sanity check plots
% for sessix = 1:length(obj)
%     for c = 1:length(obj(sessix).selectiveCells)
%         cellix = obj(sessix).selectiveCells(c);
%         psths = mySmooth(squeeze(obj(sessix).psth(:,cellix,cond2use)),31);
%         sel = mySmooth(obj(sessix).neuralselect(:,c),31);
%         subplot(1,2,1)
%         plot(obj(sessix).time,psths); legend('2AFC','AW')
%         xlim([-3 2.5])
%         subplot(1,2,2)
%         plot(obj(sessix).time,sel)
%         xlim([-3 2.5])
%         title([meta(sessix).anm meta(sessix).date '  Cell' num2str(cellix)])
%         pause
%         hold off;
%     end
% end
%% Find single-cell selectivities from reconstructed neural populations
cond2use = [12,13];
smooth = 200;

for sessix = 1:length(meta)
    space = 'null';
    rez(sessix).neuralselectivity.(space) = calcSelectivity_ReconstructedNeuralPop(sessix,obj,params,rez,space,smooth,cond2use);

    space = 'potent';
    rez(sessix).neuralselectivity.(space) = calcSelectivity_ReconstructedNeuralPop(sessix,obj,params,rez,space,smooth,cond2use);
end
%% Sanity check--for reconstructed null and potent
% cond2use = [9,10];
% for sessix = 1:length(obj)
%     for c = 1:length(obj(sessix).selectiveCells)
%         cellix = obj(sessix).selectiveCells(c);
%         % Null
%         space = 'null';
%         nullpsths = mySmooth(squeeze(rez(sessix).recon_psth.(space)(:,cellix,cond2use)),31);
%         sel = mySmooth(rez(sessix).neuralselectivity.(space)(:,c),31);
%         subplot(2,2,1)
%         plot(obj(sessix).time,nullpsths); legend('2AFC','AW')
%         xlim([-3 2.5])
%         title(space)
%         subplot(2,2,2)
%         plot(obj(sessix).time,sel)
%         xlim([-3 2.5])
%         title([meta(sessix).anm meta(sessix).date '  Cell' num2str(cellix)])
% 
%         space = 'potent';
%         nullpsths = mySmooth(squeeze(rez(sessix).recon_psth.(space)(:,cellix,cond2use)),31);
%         sel = mySmooth(rez(sessix).neuralselectivity.(space)(:,c),31);
%         subplot(2,2,3)
%         plot(obj(sessix).time,nullpsths); legend('2AFC','AW')
%         xlim([-3 2.5])
%         title(space)
%         subplot(2,2,4)
%         plot(obj(sessix).time,sel)
%         xlim([-3 2.5])
%         title([meta(sessix).anm meta(sessix).date '  Cell' num2str(cellix)])
% 
%         pause
%         hold off;
%     end
% end
%% Find correlation/selectivity explained by the reconstructed PSTHs from null or potent spaces
% Asking about the selectivty explained by the null/potent spaces only in the presample period 
for sessix = 1:length(meta)
    [nullSE,potentSE] = calcSelExplained_PSTH(sessix,obj,start,stop);
    sel_explained(sessix).null = mean(nullSE, 'omitnan');
    sel_explained(sessix).potent = mean(potentSE, 'omitnan');
end
%% Bar plot of R2 values on each session between true and predicted Null and Potent CDlate
allsel_null = [];
allsel_potent = [];
for sessix = 1:length(meta)
    allsel_null = [allsel_null, sel_explained(sessix).null];
    allsel_potent = [allsel_potent, sel_explained(sessix).potent];
end

colors = getColors_Updated();
sigcutoff = 0.05;

plotSelExplainedBar(allsel_null,allsel_potent,colors,sigcutoff)
%% Functions
function [selectivecells,spkdif] = findSelectiveCells(sessix,obj,modparams,meta,conds2comp)
currobj = obj(sessix);
nTrials = size(currobj.trialdat,3);
nCells = size(currobj.psth,2);
includedCells = [];
probenum = meta(sessix).probe;
spkstuff = currobj.clu{meta(sessix).probe};
preSampSpikes = zeros(nCells,nTrials);
for c = 1:nCells                                        % For each cell...
    cellQual = currobj.clu{probenum}(c).quality;
    % Exclude cell from analysis if it is not of proper quality
    if strcmp(cellQual,modparams.quals2excl{1}) || strcmp(cellQual,modparams.quals2excl{2}) || strcmp(cellQual,modparams.quals2excl{3})
        includedCells = [includedCells,0];
    else
        for t = 1:nTrials                                                   % Go through all of the trials
            spikeix = find(currobj.clu{probenum}(c).trial==t);              % Find the spikes for this cell that belong to the current trial
            spktms = currobj.clu{probenum}(c).trialtm_aligned(spikeix);     % Get the aligned times within the trial that the spikes occur
            prespks = length(find(spktms<samp));                            % Take the spikes which occur before the sample tone
            if ~isempty(prespks)
                preSampSpikes(c,t) = prespks;                               % Save this number
            end
        end
        includedCells = [includedCells,1];
    end
end

% Whether you are comparing the firing rates of cells across contexts or spike counts
if strcmp(modparams.measure,'FR')
    psth = obj(sessix).trialdat;                                 % Get the trial PSTH (time x cells x trials)
    presamp_psth = mean(psth(start:stop,:,:),1);                 % Take the average FR for all cells during the presamp period
    temp = squeeze(presamp_psth);                                % (cells x trials)
elseif strcmp(modparams.measure,'spkCnt')
    temp = preSampSpikes;                                        % (cells x trials)
end


for cond = 1:length(conds2comp)
    c = conds2comp(cond);
    trix = params(sessix).trialid{c};
    trix2use = randsample(trix,modparams.subTrials);
    epochAvg{cond} = temp(:,trix2use);
end
totalspksAFC = sum(epochAvg{1},2);
totalspksAW = sum(epochAvg{2},2);
spkdif = (totalspksAFC>totalspksAW);       % 1 if AFC-preferring; 0 if AW-preferring


% The p-value that you want to perform the ranksum test at
sig = 0.01;
[hyp] = getContextModulatedCells(epochAvg,sig);

selectivecells = find(hyp&includedCells);
clear epochAvg
clear includedCells
end

function selectivity = findSelectivity_PSTH(sessix,obj,del2use,cond2use,smooth)
currobj = obj(sessix);
temppsth = currobj.trialdat;
delLength = currobj.bp.ev.goCue-currobj.bp.ev.delay;
deltrix = find(delLength<(del2use+0.01)&delLength>(del2use-0.01));
psth2use = [];
for c = 1:length(cond2use)
    cond = cond2use(c);
    condtrix = params(sessix).trialid{cond};
    trix2use = deltrix(ismember(deltrix,condtrix));
    condpsth = mean(temppsth(:,:,trix2use),3,'omitnan');
    condpsth = mySmooth(condpsth,smooth);
    psth2use = cat(3,psth2use,condpsth);
end

selectivity = NaN(length(currobj.time),length(currobj.selectiveCells));
for sel = 1:length(currobj.selectiveCells)
    cellix = currobj.selectiveCells(sel);
    tempsel = psth2use(:,cellix,1)-psth2use(:,cellix,2);
    if ~currobj.spkdif(cellix)
        tempsel = -1*tempsel;
    end
    selectivity(:,sel) = mySmooth(tempsel,21);
end
end

function [nullSE,potentSE] = calcSelExplained_PSTH(sessix,obj,start,stop)
nCells = size(obj(sessix).neuralselect,2);
sel_explained(sessix).null = NaN(1,nCells);
sel_explained(sessix).potent = NaN(1,nCells);
tempnull = [];
temppot = [];
for cellix = 1:nCells
    fullsel = obj(sessix).neuralselect(start:stop,cellix);
    nullsel = rez(sessix).neuralselectivity.null(start:stop,cellix);
    potentsel = rez(sessix).neuralselectivity.potent(start:stop,cellix);

    nullR2 = corrcoef(fullsel,nullsel);
    nullR2 = nullR2(2);
    tempnull = [tempnull,nullR2];

    potentR2 = corrcoef(fullsel,potentsel);
    potentR2 = potentR2(2);
    temppot = [temppot,potentR2];

    %         subplot(1,2,1)
    %         scatter(fullsel,nullsel)
    %         xlabel('Full pop')
    %         ylabel('Null recon')
    % %         yyaxis left
    % %         plot(obj(sessix).time(start:stop),fullsel);
    % %         yyaxis right
    % %         plot(obj(sessix).time(start:stop),nullsel);
    %         title(['Null;  R2 = ' num2str(nullR2)])
    %         hold off;

    %         subplot(1,2,2)
    %         scatter(fullsel,potentsel)
    %         xlabel('Full pop')
    %         ylabel('Potent recon')
    % %         yyaxis left
    % %         plot(obj(sessix).time(start:stop),fullsel);
    % %         yyaxis right
    % %         plot(obj(sessix).time(start:stop),potentsel);
    %         title(['Potent;  R2 = ' num2str(potentR2)])
    %         hold off;
    %         pause
    nullSE = tempnull;
    potentSE = temppot;
end
end
%% Plotting functions
function plotSelExplainedBar(allsel_null,allsel_potent,colors,sigcutoff)
figure();
X = categorical({'Null','Potent'});
X = reordercats(X,{'Null','Potent'});
Y = [mean(allsel_null), mean(allsel_potent)];
b = bar(X,Y); hold on;
b.FaceColor = 'flat';
b.FaceAlpha = 0.5;
b.CData(1,:) = colors.null;
b.CData(2,:) = colors.potent;
scatter(1,allsel_null,35,[0 0 0],'filled','MarkerEdgeColor','black')
scatter(2,allsel_potent,35,[0 0 0],'filled','MarkerEdgeColor','black')
for sessix = 1:length(allsel_potent)
    plot([1,2],[allsel_null(sessix),allsel_potent(sessix)],'Color','black')
end

hyp = ttest(allsel_null,allsel_potent,'Alpha',sigcutoff);
if hyp
    scatter(1.5,1.2,30,'*','MarkerEdgeColor','black')
    title(['Presamp selectivity explained by null/pot dims; p<' num2str(sigcutoff)])
else
    title(['Presamp selectivity explained by null/pot dims; p>' num2str(sigcutoff)])
end
ylabel("R^2 value")
end