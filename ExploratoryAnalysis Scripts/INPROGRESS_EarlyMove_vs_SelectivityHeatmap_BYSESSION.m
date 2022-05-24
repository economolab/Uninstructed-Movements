% SCRIPT FOR RELATING EARLY MOVEMENTS TO SELECTIVITY ON INDIVIDUAL CELL
% BASIS

% Will generate a figure for each session:
% Panel A: Heatmap of direction selectivity for all cells recorded in the
% session (on all non-early hit trials)
% Panel B: Heatmap of direction selectivity for all cells recorded in the
% session (on only NON-EARLY MOVE hit trials)
% Panel C: Heatmap of direction selectivity for all cells recorded in the
% session (on only EARLY MOVE hit trials)

%%
clear; clc; close all;

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\ActivityModes'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Utils'));


% Saving params
outputdir = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Figures\Uninstructed Movements';
toSave = 'no';
%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % goCue or firstLick

params.lowFR               = 1; % remove clusters firing less than this val

% set conditions to use for projections
params.condition(1) = {'R&hit&~stim.enable&autowater.nums==2&~early'}; % right hits, no stim, aw off
params.condition(2) = {'L&hit&~stim.enable&autowater.nums==2&~early'}; % left hits, no stim, aw off
params.condition(3) = {'R&miss&~stim.enable&autowater.nums==2&~early'};   % error right, no stim, aw off
params.condition(4) = {'L&miss&~stim.enable&autowater.nums==2&~early'};   % error left, no stim, aw off
params.condition(5) = {'R&hit&~stim.enable&autowater.nums==1&~early'}; % right hits, no stim, aw on
params.condition(6) = {'L&hit&~stim.enable&autowater.nums==1&~early'}; % left hits, no stim, aw on
params.condition(7) = {'~hit&~miss&~stim.enable&autowater.nums==2&~early'}; % ignore, 2afc, no stim

%% SET METADATA FROM ALL RELEVANT SESSIONS/ANIMALS
meta = [];
% meta = loadJEB4_ALMVideo(meta);
% meta = loadJEB5_ALMVideo(meta);
% meta = loadJEB6_ALMVideo(meta);
meta = loadJEB7_ALMVideo(meta);
% meta = loadEKH1_ALMVideo(meta);
% meta = loadEKH3_ALMVideo(meta);
% meta = loadJGR2_ALMVideo(meta);
% meta = loadJGR3_ALMVideo(meta);

taxis = meta(end).tmin:meta(end).dt:meta(end).tmax;   % get time-axis with 0 as time of event you aligned to
taxis = taxis(1:end-1);
%% PREPROCESS DATA
objs = loadObjs(meta);

for i = 1:numel(meta)
    obj = objs{i};
    obj.condition = params.condition;
    % get trials and clusters to use
    meta(i).trialid = findTrials(obj, obj.condition);   % Get which trials pertain to the behavioral conditions you are looking at
    cluQuality = {obj.clu{meta(i).probe}(:).quality}';  % Get clusters that are of the qualities that you specified
    meta(i).cluid = findClusters(cluQuality, meta(i).quality);
    % align data
    obj = alignSpikes(obj,meta(i),params);              % Align the spike times to the event that you specified 
    % get trial avg psth, single trial data, and single trial data grouped
    % by condition (aka R 2AFC, R AW, etc.)
    obj = getPSTHs(obj,meta(i));                          
    objs{i} = obj;
end

[multi_psth,fromsess] = concatPSTH(objs);                          % Concatenate the trial-averaged PSTHs from all sessions into one structure 
%% ACTIVITY MODES

for gg = 1:length(meta)         % For all loaded sessions...
    f = figure(gg);
    f.WindowState = 'maximized';
    obj = objs{gg};
    met = meta(gg);
    smooth = 61;

    anm = obj.pth.anm;                  % Animal name
    date = obj.pth.dt;                  % Session date
    probenum = string(met.probe);       % Which probe was used  
    sesstitle = strcat(anm,{'   '},date,' ;  ','Probe ',{'   '},probenum,' ;  ','Direction Selectivity ');  % Name/title for session

    % Identify early movement trials (tongue and jaw movements during delay) 
    obj = findEarlyMoveTrials(obj);

    % Remove early movement trials from PSTH calculation
    obj = removeEarlyTrials_PSTH(obj,met);

    % PANEL A: Heat-map of R/L selectivity for all cells; all trials
    subplot(1,3,1)
    numCells = size(obj.psth,2);

    selectNorm = findDirectionSelectivity(obj, 'obj.psth',smooth);          % Find direction selectivity
    lateDelay = find(taxis>=-0.25, 1, 'first');      % t = 0.25 s before alignEvent

    sort_by = lateDelay;          % Which event/epoch you want to sort the heatmap by (choose from defined times above)
    [selectNorm,sIx] = sortSelectivity(selectNorm,sort_by,'first',[]);
    
    plotSelectivityHeatmap(taxis,sort_by,numCells,selectNorm);       % Plot heatmap for direction selectivity (for all cells in the session)
    ylabel('Neuron #','FontSize',18)
    xlabel('Time since go cue (s)','FontSize',18)
    ax = gca;
    ax.FontSize = 14;
    title('All trials','FontSize',18)
    hold off;

    % PANEL B: Heat-map of R/L selectivity for all cells --> EARLY trials
    % removed
    subplot(1,3,2)
    numCells = size(obj.psth,2);
    
    selectNorm = findDirectionSelectivity(obj, 'obj.trialpsth_noEarly',smooth);          % Find direction selectivity
    [selectNorm,sIx] = sortSelectivity(selectNorm,sort_by,'second',sIx);
    sort_by = lateDelay;          % Which event/epoch you want to sort the heatmap by (choose from defined times above)

    plotSelectivityHeatmap(taxis,sort_by,numCells,selectNorm)       % Plot heatmap for direction selectivity (for all cells in the session)
    ylabel('Neuron #','FontSize',18)
    xlabel('Time since go cue (s)','FontSize',18)
    ax = gca;
    ax.FontSize = 14;
    title('Early trials REMOVED','FontSize',18)
    hold off;

    % PANEL c: Heat-map of R/L selectivity for all cells --> EARLY trials
    % ONLY
    subplot(1,3,3)
    numCells = size(obj.psth,2);
    
    selectNorm = findDirectionSelectivity(obj, 'obj.trialpsth_Early',smooth);          % Find direction selectivity
    [selectNorm,sIx] = sortSelectivity(selectNorm,sort_by,'second',sIx);
    sort_by = lateDelay;          % Which event/epoch you want to sort the heatmap by (choose from defined times above)

    plotSelectivityHeatmap(taxis,sort_by,numCells,selectNorm)       % Plot heatmap for direction selectivity (for all cells in the session)
    ylabel('Neuron #','FontSize',18)
    xlabel('Time since go cue (s)','FontSize',18)
    ax = gca;
    ax.FontSize = 14;
    title('Early trials ONLY','FontSize',18)
    hold off;

    sgtitle(sesstitle)

     % Save the figure to the output directory
    if strcmp(toSave,'yes')
        saveas(gcf,fullfile(outputdir,sesstitle),'jpeg')
        close all
    end
   
end