% Finding alignment of single cells to CDContext in the null or potent subspaces from neural activity that resides within the Null and Potent spaces
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--- CDContext found usual way (all DR - WC trials used and averaged together) ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

load([basepth '\Uninstructed-Movements\ContextColormap.mat']);
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
meta = loadJEB19_ALMVideo(meta,datapth);

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
% -----------------------------------------------------------------------

for sessix = 1:numel(meta)

    % -- input data
     trialdat_zscored = zscore_singleTrialNeuralData(obj(sessix));
    zscored(sessix).trialdat =  trialdat_zscored;

    % -- Calculate the null and potent spaces for each session
    cond2use = [2 3 4 5];   % All 2AFC hit trials, all AW hit trials (NUMBERING ACCORDING TO PARAMS.CONDITION)
    nullalltime = 0;        % use all time points to estimate null space if 1
    AWonly = 0;             % use only AW to find null and potent spaces 
    delayOnly = 0;          % use only delay period to find null and potent spaces
    cond2proj = [2 3];       % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    rez(sessix) = singleTrial_elsayed_np(trialdat_zscored, obj(sessix), me(sessix), params(sessix), cond2use, cond2proj,nullalltime,AWonly,delayOnly);

    % -- Find coding dimensions from RECONSTRUCTED full neural activity which is reconstructed from the null and potent spaces
    cond2use = [1 2];            % (NUMBERING ACCORDING TO THE CONDITIONS PROJECTED INTO NULL AND POTENT SPACES, i.e. which of the conditions specified in 'cond2proj' above do you want to use?)
    cond2proj = [1 2];           % 2AFC hits, AW hits, 2AFC miss, AW miss (corresponding to null/potent psths in rez)
    cond2use_trialdat = [2 3];   % (NUMBERING ACCORDING TO PARAMS.CONDITION)
    cd_null(sessix) = getCodingDimensions_Context(rez(sessix).recon_psth.null,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);
    cd_potent(sessix) = getCodingDimensions_Context(rez(sessix).recon_psth.potent,trialdat_zscored,obj(sessix),params(sessix),cond2use,cond2use_trialdat, cond2proj);

end
%% Project single trials onto Null and Potent CDs
disp('----Projecting single trials onto CDContext----')
cd = 'context';

[cd_null,cd_potent] = getNPSingleTrialProjs(obj,cd,cd_null,cd_potent,rez); 
%% find DR/WC selective cells per session
% only using cells with significant selectivity in this analysis
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
edges = [trialstart samp];
cond2use = [6 7];
for i = 1:numel(obj)
    cluix{i} = findSelectiveCells(obj(i),params(i),edges,cond2use);
end

%% reconstruct single cell activity from CDcontext projs in either null or potent space

cdix = 1; % index of cdchoice in the projs
cond2use = {'hit|miss'}; % specifying which trials 
for sessix = 1:numel(meta)
    clear trialdat W proj 
    disp(['Session ' num2str(sessix) '/' num2str(numel(meta))])
    trix = findTrials(obj(sessix), cond2use);
    trix = trix{1};
    clus = find(cluix{sessix});

    % get full single trial data (will compare reconstructed against this)
%     trialdat.full = zscore_singleTrialNeuralData(obj(sessix));
    trialdat.full = permute(obj(sessix).trialdat(:,:,trix),[1 3 2]);
%     trialdat.full = trialdat.full(:,trix,:); % (time,trials,neurons);

    % get CDs
    W.null = cd_null(sessix).cd_mode_orth(:,cdix);
    W.potent = cd_potent(sessix).cd_mode_orth(:,cdix);
    fns = {'null','potent'};
    for j = 1:numel(fns)
        % single trials neural activity reconstructed from n/p
        trialdat.(fns{j}) = rez(sessix).recon.(fns{j})(:,trix,:); % (time,trials,dims)
        % project onto Wcontext
        proj.(fns{j}) = tensorprod(trialdat.(fns{j}),W.(fns{j}),3,1);
        % reconstruct data from CD context proj
        trialdat.recon.(fns{j}) = tensorprod(proj.(fns{j}),W.(fns{j}),3,2);

        % for each cell, get R^2 b/w it's original data and reconstructed
        for k = 1:numel(clus) % for each cell
            thisclu = clus(k);
            % calculcate variance explained by CD context
            orig = trialdat.full(:,:,thisclu); % (time,trials)

            fr = mean(mean(orig)); % subspace contribution method
            % weight = norm(W.(fns{j})(k));
            % r2.(fns{j}){sessix}(k) = fr*weight;

            recon = trialdat.recon.(fns{j})(:,:,thisclu); % (time,trials) % ve by recon method
            mdl = fitlm(orig(:),recon(:));
            r2.(fns{j}){sessix}(k) = mdl.Rsquared.Ordinary;
          
        end
    end
end


%% plot
close all

% concatenate R^2s from null and potent spaces into two vectors
alln = [];
allp = [];
for sessix = 1:numel(meta)
    n = r2.null{sessix};
    alln = [alln n];
    p = r2.potent{sessix};
    allp = [allp p];
    % scatter(n,p,10,'filled','MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5]);
    
    trix = findTrials(obj(sessix), cond2use);
    trix = trix{1};
    trialdat.full = permute(obj(sessix).trialdat(:,:,trix),[1 3 2]);

    clus = find(cluix{sessix});
    %%%%%%%%---------------SANITY CHECK------------------%%%%%%%%%%%%%%%%
    for k = 1:numel(clus) % for each cell
        thisclu = clus(k);
        % calculcate variance explained by CD context
        
        orig = trialdat.full(:,:,thisclu); % (time,trials)

%         tempAL = (n(k) - p(k)) ./ (p(k) + n(k)); % calculate alignment index
%         subplot(1,2,1)
%         imagesc(orig');
%         ylabel('Trials')
%         subplot(1,2,2)
%         plot(obj(sessix).time, mySmooth(obj(sessix).psth(:,thisclu,6),31)); hold on; plot(obj(sessix).time,mySmooth(obj(sessix).psth(:,thisclu,7),31)); hold off; 
%         xline(0,'k--')
%         xlabel('Time from go cue')
%         ylabel('FR')
%         legend('DR','WC')
%         sgtitle(num2str(tempAL))
%         pause
    %%%%%%%%%%-------------------------------------------%%%%%%%%%%%%%%%%%%%
    end
end
alignment = (alln - allp) ./ (allp + alln); % calculate alignment index

% histogram
% f = figure;
% f.Position = [644   483   338   231];
% ax = gca;
% f.Renderer = 'painters';
% ax = prettifyPlot(ax);
hold on;

h = histogram(alignment,40,'edgecolor','none','Normalization','count');
ylabel('# Neurons')
xlabel('CDContext alignment')
xline(0,'k--')





