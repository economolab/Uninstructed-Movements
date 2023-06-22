% Big context-selectivity heatmap of single-neurons

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
addpath(genpath(fullfile(utilspth,'fig3')));
figpth = [basepth '\Uninstructed-Movements\Fig 3'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context switching')));
figpth = [basepth '\Uninstructed-Movements\Fig 6'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Context_funcs')));
figpth = [basepth '\Uninstructed-Movements\Fig 5'];
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
params.condition(1) = {'hit&~early&~stim.enable&~autowater'};                          %  hit, no stim, aw off
params.condition(end+1) = {'hit&~early&~stim.enable&autowater'};                          %  hit, no stim, aw on

params.tmin = -3;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 30;

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
%% Concat PSTHs from all cells across sessions
psth_all = [];
for sessix = 1:length(meta)
    temp = obj(sessix).psth;
    psth_all = cat(2,psth_all, temp);
end
% Find context selectivity (2AFC - AW)
selectivity = squeeze(psth_all(:,:,1)-psth_all(:,:,2));     % Avg PSTH on 2AFC trials - Avg PSTH on AW trials

% Do normalization (as in Li et al, 2015 'A motor cortex circuit...')
maxSelect =  max(abs(selectivity),[],1);                    % Max magnitude of context selectivity for each cell
selectNorm = NaN(length(obj(1).time),length(maxSelect));
for c = 1:length(maxSelect)                                 % For every cell... 
    selectNorm(:,c) = selectivity(:,c)./maxSelect(c);       % Normalize all selectivity values to the max selectivity val
end                                                         % Each cell will have selectivity values -1 to 1
%% Sort by pre-sample selectivity
% Find the times corresponding to trial start and the sample period
trialstart = median(obj(1).bp.ev.bitStart)-median(obj(1).bp.ev.(params(1).alignEvent));
start = find(obj(1).time>trialstart,1,'first');
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time<samp,1,'last');

% Sort the selectivity PSTHs
temp = mean(selectNorm(start:stop,:),1,'omitnan');          % Find the average presample selectivity for each neuron
[~,sortix] = sort(temp,'descend');                          % Sort the selectivity PSTHs in descending order according to presamp selectivity
sorted_select = selectNorm(:,sortix);
%% Find cells that are significantly modulated by context (ranksum test, p-value = 0.01)
% As in Inagaki et al., Cell, 2022 ('A midbrain-thalamus...')
% modulatedCells = (1 x nCells) array where 1 means the cell is context-selective and 0 means it is not
modulatedCells = [];
modparams.quals2excl = {'Poor','Multi','Noisy','Garbage'};
includedCells = [];
modparams.sm = 30;                                                    % Amount that you want to smooth PSTHs by
modparams.measure = 'FR';                                             % Whether you want to compare firing rate ('FR') or spike counts ('spkCnt') 
modparams.subTrials = 35;
for sessix = 1:length(meta)                                 % For every session...   
    currobj = obj(sessix);
    nTrials = size(currobj.trialdat,3);
    nCells = size(currobj.psth,2);
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
    
    for cond = 1:size(obj(sessix).psth,3)
        trix = params(sessix).trialid{cond};
        trix2use = randsample(trix,modparams.subTrials);
        epochAvg{cond} = temp(:,trix2use);
    end
    

    % The p-value that you want to perform the ranksum test at
    sig = 0.01;
    [hyp] = getContextModulatedCells(epochAvg,sig);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for c = 1:nCells                                        % For each cell...
%         cellQual = currobj.clu{probenum}(c).quality;
%         if ~strcmp(cellQual,quals2excl{1}) || ~strcmp(cellQual,quals2excl{2}) || ~strcmp(cellQual,quals2excl{3})
%             subplot(1,2,1)
%             nBins = 20;
%             histogram(epochAvg{1}(c,:),nBins); hold on; 
%             histogram(epochAvg{2}(c,:),nBins); hold off;
%             title(num2str(hyp(c)))
% 
%             subplot(1,2,2)
%             plot(obj(sessix).psth(:,c,1)); hold on; plot(obj(sessix).psth(:,c,2));
%             hold off;
%             pause
%         end
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modulatedCells = [modulatedCells, hyp];


    clear epochAvg
end
%% 
ctxtSelective = sum(includedCells&modulatedCells);
nIncl = sum(includedCells);
pctSelective = 100*(ctxtSelective/nIncl);

disp('---Summary Statistics for number of context selective neurons---')
disp([num2str(pctSelective) ' % of cells show context-selectivity in the presample period (' num2str(ctxtSelective) '/' num2str(nIncl) ' cells)'])
t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
disp(t)
%% Sort according to how we want to plot it for figure
clearvars -except obj meta params includedCells modulatedCells selectNorm psth_all start stop 

% Get the indices of all cells which are included and context-modulated
tempix = find(includedCells&modulatedCells);
selectivityMod = selectNorm(:,tempix);                                  % Selectivity values for these cells
psthMod = psth_all(:,tempix,:);                                         % PSTHs for these cells                
clear tempix 

% Get the presample average selectivity for all modulated cells
preAvg = mean(selectivityMod(start:stop,:),1,'omitnan');

% Find cells which are selective for 2AFC
afcix = find(preAvg>0);                                                 % Cells with a positive selectivity value
included.selectivity.AFC = selectivityMod(:,afcix);
included.psth.AFC = psthMod(:,afcix,:);
temp = mean(included.selectivity.AFC(start:stop,:),1,'omitnan');        % Get the average pre-sample selectivity for 2AFC selective cells
[~,sortix] = sort(temp,'descend');                                      % Sort 2AFC-selective cells in order of pre-sample selectivity
included.selectivity.AFC = included.selectivity.AFC(:,sortix);
included.psth.AFC = included.psth.AFC(:,sortix,:);                      % Sort PSTHs in same order

% Do the same thing for AW-preferring cells
awix = find(preAvg<0);                                                  % Cells with a negative selectivity value are AW-selective
included.selectivity.AW = selectivityMod(:,awix);
included.psth.AW = psthMod(:,awix,:);
temp = mean(included.selectivity.AW(start:stop,:),1,'omitnan');
[~,sortix] = sort(temp,'ascend');  
included.selectivity.AW = included.selectivity.AW(:,sortix);
included.psth.AW = included.psth.AW(:,sortix,:);

% Concatenate the sorted selectivity and PSTH values for 2AFC and AW preferring cells (preserving the sorted order)
selectivityMod = [included.selectivity.AFC,included.selectivity.AW];
psthMod = cat(2,included.psth.AFC,included.psth.AW);

% Number of cells in each type
nums.AFC = size(included.selectivity.AFC,2);
nums.AW = size(included.selectivity.AW,2);

% Get the cells which are not modulated by context and sort them according to presample selectivity
tempix = find(includedCells&~modulatedCells);
selectivityNonMod = selectNorm(:,tempix);
psthNonMod = psth_all(:,tempix,:);
% Sort non-modulated cells by presample selectivity
% temp = mean(selectivityNonMod(start:stop,:),1,'omitnan');
% [~,sortix] = sort(temp,'descend');  
% selectivityNonMod = selectivityNonMod(:,sortix);
% psthNonMod = psthNonMod(:,sortix,:);
% Randomize order of non-modulated cells
randix = randsample(1:size(selectivityNonMod,2),size(selectivityNonMod,2));
selectivityNonMod = selectivityNonMod(:,randix);
psthNonMod = psthNonMod(:,randix,:);

% Concatenate such that cells are ordered with: 2AFC-selective first (sorted), AW-selective, then non-selective
toplot_all = [selectivityMod,selectivityNonMod];
psthplot_all = cat(2,psthMod, psthNonMod);
%%
disp('---Summary Statistics for types of context selective neurons---')
disp([num2str(nums.AFC) ' cells out of ' num2str(sum([nums.AFC,nums.AW])) ' are 2AFC-preferring'])
disp([num2str(nums.AW) ' cells out of ' num2str(sum([nums.AFC,nums.AW])) ' are AW-preferring'])
t = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss Z');
disp(t)
%%
% for c = 1:size(psthplot_all,2)                                 % For every cell... 
%     temp = squeeze(psthplot_all(:,c,:));                                % Get its max firing rate across both conditions
%     maxFR = max(abs(temp));
%     psthplot_norm(:,c,:) = temp./maxFR;                 % Normalize all FRs to the max FR for that cell
% end 
%%
clearvars -except obj meta params includedCells modulatedCells toplot_all psthplot_all nums selectNorm psth_all start stop

load('C:\Code\Uninstructed-Movements\ContextColormap.mat');

tLim.start = -3;
tLim.stop = 1;
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.(params(1).alignEvent));
nCells = size(toplot_all,2);

subplot(1,3,1)
imagesc(obj(1).time,1:nCells,psthplot_all(:,:,1)'); hold on
colorbar
ax = gca;
colormap(ax,flipud(ContextColormap))
clim([0 70])
line([tLim.start,tLim.stop],[nums.AFC,nums.AFC],'Color','black','LineWidth',1.25)
line([tLim.start,tLim.stop],[nums.AFC+nums.AW,nums.AFC+nums.AW],'Color','black','LineWidth',1.25)
line([samp,samp],[1,nCells],'Color','black','LineStyle','--','LineWidth',1.5)
xlim([tLim.start tLim.stop])
title('2AFC Trials')
ylabel('Neuron Number')

subplot(1,3,2)
imagesc(obj(1).time,1:nCells,psthplot_all(:,:,2)'); hold on
colorbar
ax = gca;
colormap(ax,flipud(ContextColormap))
clim([0 70])
line([tLim.start,tLim.stop],[nums.AFC,nums.AFC],'Color','black','LineWidth',1.25)
line([tLim.start,tLim.stop],[nums.AFC+nums.AW,nums.AFC+nums.AW],'Color','black','LineWidth',1.25)
line([samp,samp],[1,nCells],'Color','black','LineStyle','--','LineWidth',1.5)
xlim([tLim.start tLim.stop])
xlabel('Time from go cue/water drop (s)')
title('AW Trials')

subplot(1,3,3)
imagesc(obj(1).time,1:nCells,toplot_all'); hold on
colorbar
ax = gca;
colormap(ax,flipud(ContextColormap))
line([tLim.start,tLim.stop],[nums.AFC,nums.AFC],'Color','black','LineWidth',1.25)
line([tLim.start,tLim.stop],[nums.AFC+nums.AW,nums.AFC+nums.AW],'Color','black','LineWidth',1.25)
line([samp,samp],[1,nCells],'Color','black','LineStyle','--','LineWidth',1.5)
xlim([tLim.start tLim.stop])
title('Selectivity')
%% Plot all cells
% nCells = size(sorted_select,2);
% imagesc(obj(1).time,1:nCells,sorted_select')
% colorbar
% colormap("jet")
% %clim([-0.2 0.2])
% hold on;
% line([samp,samp],[1,nCells],'Color','black','LineStyle','--','LineWidth',1.5)
% xlabel('Time from water drop (s)')
% ylabel('Neuron Number')
%%
function [hyp] = getContextModulatedCells(epochAvg,sig)
% Produces a (2 x nCells) array of p-values and hyp test results for every
% cell in L and R trials
                                                
    nCells = size(epochAvg{1},1);                   % Get the number of cells
    pvals = zeros(1,nCells);                     % Store p-values for each cell
    hyp = zeros(1,nCells);                       % Store hyp test results for each cell
    for c = 1:nCells
        [pR,hR] = ranksum(epochAvg{1}(c,:),epochAvg{2}(c,:),'alpha',sig);     % Run Wilcoxon Ranksum test on each cell between 2AFC and AW conditions
        pvals(1,c) = pR;                                                % Store p-values
        hyp(1,c) = hR;                                                  % Store outcome of hypothesis test
    end 

end   % getContextModCells