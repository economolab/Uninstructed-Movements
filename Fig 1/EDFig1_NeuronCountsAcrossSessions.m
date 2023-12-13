% ED Figure 1 -- Counts of single and multi-units across sessions
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
addpath(genpath(fullfile(utilspth,'fig1')));
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
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off
params.condition(end+1) = {'R&no&~stim.enable&~autowater&~early'};              % no right, no stim, aw off
params.condition(end+1) = {'L&no&~stim.enable&~autowater&~early'};              % no left, no stim, aw off
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};               % all hits, no stim, aw off


params.tmin = -2.5;
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

meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written

%% LOAD DATA

% ----------------------------------------------
% -- Neural Data --
% obj (struct array) - one entry per session
% params (struct array) - one entry per session
% ----------------------------------------------
[obj,params] = loadSessionData(meta,params);
%% Get number of single neurons and multi units for each session
numSingle = NaN(1,length(meta));            % [1 x nSessions]
numMulti = NaN(1,length(meta));             % [1 x nSessions]
for sessix = 1:length(meta)                 % For every session...
    nProbes = length(meta(sessix).probe);
    sessSingle = [];
    sessMulti = [];
    for num = 1:nProbes                     % For each probe...
        probenum = meta(sessix).probe(num);         % Get the probe number
        clus = obj(sessix).clu{probenum};           % Get the cluster information from that probe
        singlemask = zeros(1,length(clus));         % [1 x nCells]
        multimask = zeros(1,length(clus));          % [1 x nCells]
        for cellix = 1:length(clus)                   % For each cell...
            cellqual = clus(cellix).quality;              % Get the quality of that cell
            if contains(cellqual,'air')|| contains(cellqual,'ood') ...          % If fair, good, excellent or great...
                    || contains(cellqual,'xcellent')|| contains(cellqual,'reat')
                singlemask(cellix) = 1;                       % Count it as a single unit
            elseif contains(cellqual,'ulti')            % If multi...
                multimask(cellix) = 1;                        % Count it as a multi unit
            end
        end
        sessSingle = [sessSingle, sum(singlemask)]; % Concatenate across probes if needed
        sessMulti = [sessMulti, sum(multimask)];
    end
    if size(sessSingle,2)==1                % If one probe...
        numSingle(sessix) = sessSingle;         % Number of single neurons for this session 
        numMulti(sessix) = sessMulti;
    elseif size(sessSingle,2)==2            % If two probes...
        numSingle(sessix) = sum(sessSingle);    % The num of single neurons for the session is the sum across probes
        numMulti(sessix) = sum(sessMulti);
    end
end
%% Make stacked bar plot to illustrate num neurons per session
% Height of each bar = sum of elements in row

sessnumDRWC = 12;

xticklabs = cell(1,length(meta));
for sessix = 1:length(meta)
    label = [meta(sessix).anm ';' meta(sessix).date];
    xticklabs{sessix} = label;
end

y = [numSingle', numMulti'];
bar(y,'stacked')
xline(sessnumDRWC+0.515,'k--')
xticks(1:length(meta))
xticklabels(xticklabs)
xtickangle(45)
set(gca,'TickDir','out')
legend('Single units','Multi units','Location','best')
ylabel('Number of neurons')
xlabel('Session #')
title('Fixed delay sessions')