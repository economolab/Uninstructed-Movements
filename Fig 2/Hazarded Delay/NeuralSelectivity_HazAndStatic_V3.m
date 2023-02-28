% From Inagaki 2019 Nature (Discrete Attractor Dynamics paper)
% "Neurons with significant selectivity (two-sided Wilcoxon rank sum test comparing spike counts in two correct trial types, P < 0.05) during the delay epochs were classified as selective cells. Selective cells were classified into lick-right preferring versus lick-left preferring, on the basis of their total spike counts during the delay epoch.
% For peri-stimulus time histograms (Figs. 5, 6, Extended Data Fig. 8k), only correct trials were included. For the peri-stimulus time histograms and selectivity of the 
% random delay task (Fig. 6, Extended Data Fig. 8k), only spikes before the go cue were pooled. Spikes were averaged over 100 ms with a 1-ms sliding window."

clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Munib Uninstruct Move\uninstructedMovements_v2';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
figpth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\Fig 3';
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Hazarded Delay')));
%% PARAMETERS
params.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for                            % all trials
params.condition(1) = {'R&hit&~stim.enable&~autowater&~early'};             % R 2AFC hits, no stim
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % L 2AFC hits, no stim

params.tmin = -1.6;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 31;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% vid features to use
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)
params.advance_movement = 0;

% Haz delay params
params.delay(1) = 0.3000;
params.delay(2) = 0.6000;
params.delay(3) = 1.2000;
params.delay(4) = 1.8000;
params.delay(5) = 2.4000;
params.delay(6) = 3.6000;
%% SPECIFY DATA TO LOAD

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

meta = [];

% --- ALM --- 
meta = loadJEB11_ALMVideo(meta,datapth);
meta = loadJEB12_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written
%% LOAD DATA
% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);
%% Find which cells are selective for trial-type during the delay period 
%clearvars -except obj meta params
smooth = 50;
modparams.quals2excl = {'Poor','Noisy'};
modparams.sm = 30;                                                    % Amount that you want to smooth PSTHs by
modparams.subTrials = 40;
for sessix = 1:length(meta)
    includedCells = [];
    currobj = obj(sessix);
    nTrials = size(currobj.trialdat,3);
    nCells = size(currobj.psth,2);
    probenum = meta(sessix).probe;
    DelaySpikes = zeros(nCells,nTrials);
    for c = 1:nCells                                        % For each cell...
        cellQual = currobj.clu{probenum}(c).quality;
        % Exclude cell from analysis if it is not of proper quality 
        if strcmp(cellQual,modparams.quals2excl{1}) || strcmp(cellQual,modparams.quals2excl{2})
            includedCells = [includedCells,0];
        else
            for t = 1:nTrials                                                   % Go through all of the trials
                spikeix = find(currobj.clu{probenum}(c).trial==t);              % Find the spikes for this cell that belong to the current trial
                spktms = currobj.clu{probenum}(c).trialtm_aligned(spikeix);     % Get the aligned times within the trial that the spikes occur
                del = currobj.bp.ev.delay(t)-currobj.bp.ev.(params(sessix).alignEvent)(t);
                go = currobj.bp.ev.goCue(t)-currobj.bp.ev.(params(sessix).alignEvent)(t);
                delspks = length(find(spktms<go&spktms>del));                            % Take the spikes which occur before the sample tone
                if ~isempty(delspks)
                    DelaySpikes(c,t) = delspks;                                 % Save this number
                end
            end
            includedCells = [includedCells,1];
        end
    end

    temp = DelaySpikes;                                        % (cells x trials)
    
    % Sub-sample trials
    for cond = 1:size(obj(sessix).psth,3)
        trix = params(sessix).trialid{cond};
        trix2use = randsample(trix,modparams.subTrials);
        epochAvg{cond} = temp(:,trix2use);
    end
    totalspksR = sum(epochAvg{1},2);
    totalspksL = sum(epochAvg{2},2);
    spkdif = (totalspksR>totalspksL);       % 1 if R-preferring; 0 if L-preferring
   
    % The p-value that you want to perform the ranksum test at
    sig = 0.05;
    [selectiveCells] = getSelectiveCells(epochAvg,sig);

    obj(sessix).selectiveCells = find(includedCells&selectiveCells);
    obj(sessix).spkdif = spkdif; 
end
%% Calculate selectivity (spike rate difference) for all selective cells (on a given delay duration)
del2use = 1.2000;
cond2use = [1,2];
smooth = 200;
selectivity_All = [];
for sessix = 1:length(meta)
    currobj = obj(sessix);
    temppsth = currobj.trialdat;
    delLength = currobj.bp.ev.goCue-currobj.bp.ev.delay;
    deltrix = find(delLength<(del2use+0.01)&delLength>(del2use-0.01));
    psth2use = [];
    for c = cond2use
        condtrix = params(sessix).trialid{c};
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
    selectivity_All = [selectivity_All,selectivity];
end
%%
tempobj = obj(1); tempparams = params(1);
clearvars -except selectivity_All tempobj tempparams
obj = tempobj; params = tempparams;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------Static Delay-----------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PARAMETERS
ctrlparams.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
ctrlparams.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
ctrlparams.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

ctrlparams.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for                            % all trials
ctrlparams.condition(1) = {'R&hit&~stim.enable&~autowater&~early'};             % R 2AFC hits, no stim
ctrlparams.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % L 2AFC hits, no stim

ctrlparams.tmin = -1.6;
ctrlparams.tmax = 2.5;
ctrlparams.dt = 1/200;

% smooth with causal gaussian kernel
ctrlparams.smooth = 31;

% cluster qualities to use
ctrlparams.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% vid features to use
ctrlparams.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

ctrlparams.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance
ctrlparams.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)
ctrlparams.advance_movement = 0;

datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';

ctrlmeta = [];

% --- ALM --- 
ctrlmeta = loadJEB6_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB7_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadEKH1_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadEKH3_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJGR2_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJGR3_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB14_ALMVideo(ctrlmeta,datapth);
ctrlmeta = loadJEB15_ALMVideo(ctrlmeta,datapth);

ctrlparams.probe = {ctrlmeta.probe}; % put probe numbers into params, one entry for element in ctrlmeta, just so i don't have to change code i've already written
%% LOAD DATA
% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[ctrlobj,ctrlparams] = loadSessionData(ctrlmeta,ctrlparams);
%% Find which cells are selective for trial-type during the delay period 
smooth = 50;
modparams.quals2excl = {'Poor','Noisy'};
modparams.sm = 30;                                                    % Amount that you want to smooth PSTHs by
modparams.subTrials = 40;
for sessix = 1:length(ctrlmeta)
    includedCells = [];
    currobj = ctrlobj(sessix);
    nTrials = size(currobj.trialdat,3);
    nCells = size(currobj.psth,2);
    probenum = ctrlmeta(sessix).probe;
    if length(probenum)>1
        probenum = 3;
        nCells1 = length(currobj.clu{1});
        nCells2 = length(currobj.clu{2});
        currobj.clu{3} = currobj.clu{1};
        for i = 1:nCells2
            currobj.clu{3}(nCells1+i).quality = currobj.clu{2}(i).quality;
            currobj.clu{3}(nCells1+i).trial = currobj.clu{2}(i).trial;
            if ~isfield(currobj.clu{2}(1),'trialtm_aligned')
                for clu = 1:numel(currobj.clu{2})
                    event = currobj.bp.ev.(ctrlparams(sessix).alignEvent)(currobj.clu{2}(clu).trial);
                    currobj.clu{2}(clu).trialtm_aligned = currobj.clu{2}(clu).trialtm - event;
                end
            end
            currobj.clu{3}(nCells1+i).trialtm_aligned = currobj.clu{2}(i).trialtm_aligned;
        end
        
    end

    DelaySpikes = zeros(nCells,nTrials);
    for c = 1:nCells                                        % For each cell...
        cellQual = currobj.clu{probenum}(c).quality;
        % Exclude cell from analysis if it is not of proper quality 
        if strcmp(cellQual,modparams.quals2excl{1}) || strcmp(cellQual,modparams.quals2excl{2})
            includedCells = [includedCells,0];
        else
            for t = 1:nTrials                                                   % Go through all of the trials
                spikeix = find(currobj.clu{probenum}(c).trial==t);              % Find the spikes for this cell that belong to the current trial
                if ~isempty(spikeix)
                    spktms = currobj.clu{probenum}(c).trialtm_aligned(spikeix);     % Get the aligned times within the trial that the spikes occur
                    del = currobj.bp.ev.delay(t)-currobj.bp.ev.(ctrlparams(sessix).alignEvent)(t);
                    go = currobj.bp.ev.goCue(t)-currobj.bp.ev.(ctrlparams(sessix).alignEvent)(t);
                    delspks = length(find(spktms<go&spktms>del));                            % Take the spikes which occur before the sample tone
                    if ~isempty(delspks)
                        DelaySpikes(c,t) = delspks;                                 % Save this number
                    end
                end
            end
            includedCells = [includedCells,1];
        end
    end

    temp = DelaySpikes;                                        % (cells x trials)
    
    % Sub-sample trials
    for cond = 1:size(ctrlobj(sessix).psth,3)
        trix = ctrlparams(sessix).trialid{cond};
        trix2use = randsample(trix,modparams.subTrials);
        epochAvg{cond} = temp(:,trix2use);
    end
    totalspksR = sum(epochAvg{1},2);
    totalspksL = sum(epochAvg{2},2);
    spkdif = (totalspksR>totalspksL);       % 1 if R-preferring; 0 if L-preferring
   
    % The p-value that you want to perform the ranksum test at
    sig = 0.05;
    [selectiveCells] = getSelectiveCells(epochAvg,sig);

    ctrlobj(sessix).selectiveCells = find(includedCells&selectiveCells);
    ctrlobj(sessix).spkdif = spkdif; 
end
%% Calculate selectivity (spike rate difference) for all selective cells (on a given delay duration)
del2use = 0.9;
cond2use = [1,2];
smooth = 200;
selectivity_AllCtrl = [];
for sessix = 1:length(ctrlmeta)
    currobj = ctrlobj(sessix);
    temppsth = currobj.trialdat;
    delLength = currobj.bp.ev.goCue-currobj.bp.ev.delay;
    deltrix = find(delLength<(del2use+0.01)&delLength>(del2use-0.01));
    psth2use = [];
    for c = cond2use
        condtrix = ctrlparams(sessix).trialid{c};
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
    selectivity_AllCtrl = [selectivity_AllCtrl,selectivity];
end
%%
clearvars -except selectivity_All selectivity_AllCtrl obj meta params kin ctrlobj ctrlmeta ctrlparams ctrlkin

colors = getColors_Updated();
ctrlcol = colors.afc;
hazcol = [0.5 0.5 0.5];
go = mode(ctrlobj(1).bp.ev.goCue)-mode(ctrlobj(1).bp.ev.(ctrlparams(1).alignEvent));
ctrlstop = find(ctrlobj(1).time<go,1,'last');
del2use = 1.2;
smooth = 70;
alph = 0.2;

figure();
toplot = mean(selectivity_All,2,'omitnan');
%err = 1.96*(std(selectivity_All,0,2,'omitnan')/sqrt(size(selectivity_All,2)));
err = std(selectivity_All,0,2,'omitnan')/sqrt(size(selectivity_All,2));
ax = gca;
shadedErrorBar(obj(1).time, toplot, err ,{'Color',hazcol,'LineWidth',2.5}, alph, ax); hold on;

toplot = mean(selectivity_AllCtrl,2,'omitnan');
%err = 1.96*(std(selectivity_AllCtrl,0,2,'omitnan')/sqrt(size(selectivity_AllCtrl,2)));
err = std(selectivity_AllCtrl,0,2,'omitnan')/sqrt(size(selectivity_AllCtrl,2));
ax = gca;
shadedErrorBar(obj(1).time(1:ctrlstop), toplot(1:ctrlstop), err(1:ctrlstop) ,{'Color',ctrlcol,'LineWidth',2.5}, alph, ax);
%plot(obj(1).time,toplot,'Color',ctrlcol,'LineStyle','--','LineWidth',1.75)

xline(0,'LineStyle','--','Color','black','LineWidth',1.5)
xline(go,'LineStyle','--','Color',ctrlcol,'LineWidth',1.5)
xline(del2use,'LineStyle','--','Color',hazcol,'LineWidth',1.5)
%xline(-1.3,'LineStyle','-.','Color',[0.75 0.75 0.75],'LineWidth',1.5)
%legend(' ',' ',' ','Haz',' ',' ',' ','Static','Location','best')
xlim([-1.3 del2use])
xlabel('Time from delay onset (s)')
ylabel('Selectivity (spks/s)')