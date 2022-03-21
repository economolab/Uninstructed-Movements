% SCRIPT FOR RELATING CHOICE MODE ACTIVITY TO JAW VELOCITY IN A GRADED
% MANNER (HOW DOES CHOICE MODE CHANGE WHEN JAW VELOCITY CHANGES?)

% Will generate summary figure relating jaw velocity to choice coding
% direction

% PARAMS THAT NEED TO BE SET:
% Saving params: 'toSave' --> whether or not to save figures
% Run params: 'earlytrials' --> which method you want to use to identify
% early movement trials
% 'moveThresh' --> what percentage of the delay period you want to be used
% for identifying early move trials
% 'alignEvent' --> which behavioral event you want to align the PSTHs and
% modes to
%%
clear; clc; close all;

addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\ActivityModes'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Data-Loading-Scripts'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements'));
addpath(genpath('C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Utils'));

% addpath(genpath('C:\Code\ActivityModes'));
% addpath(genpath('C:\Code\Uninstructed Movements\Uninstructed-Movements\DataLoadingScripts'));
% addpath(genpath('C:\Code\Uninstructed-Movements'));
% addpath(genpath('C:\Code\Utils'));


% Saving params
outputdir = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Figures\Uninstructed Movements';
toSave = 'yes';
%% SET RUN PARAMS

% Which method you want to use to identify early movement trials:
% 'motionEnergy' or 'DeepLabCut'
params.earlytrials         =  'motionEnergy';
params.moveThresh          = 0.15;      % What percentage of the delay period you want to use for identifying early move trials
params.alignEvent          = 'goCue';   % goCue or firstLick
params.lowFR               = 1; % remove clusters firing less than this val
params.dt = 0.05;

% set conditions to use for projections
params.condition(1) = {'R&hit&~stim.enable&autowater.nums==2&~early'}; % right hits, no stim, aw off
params.condition(2) = {'L&hit&~stim.enable&autowater.nums==2&~early'}; % left hits, no stim, aw off
params.condition(3) = {'R&miss&~stim.enable&autowater.nums==2&~early'};   % error right, no stim, aw off
params.condition(4) = {'L&miss&~stim.enable&autowater.nums==2&~early'};   % error left, no stim, aw off
params.condition(5) = {'R&hit&~stim.enable&autowater.nums==1&~early'}; % right hits, no stim, aw on
params.condition(6) = {'L&hit&~stim.enable&autowater.nums==1&~early'}; % left hits, no stim, aw on
params.condition(7) = {'~hit&~miss&~stim.enable&autowater.nums==2&~early'}; % ignore, 2afc, no stim


% set conditions used for finding the modes
aw = '2'; % 1-on, 2-off
stim = '0'; % 0-off
params.modecondition(1) = {['R&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};     % R hits, 2afc, stim on/off, not early
params.modecondition(2) = {['L&hit&autowater.nums==' aw '&stim.num==' stim '&~early']};     % L hits, 2afc, stim on/off, not early
params.modecondition(3) = {['R&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};    % R miss, 2afc, stim on/off, not early
params.modecondition(4) = {['L&miss&autowater.nums==' aw '&stim.num==' stim '&~early']};    % L miss, 2afc, stim on/off, not early
params.modecondition(5) = {['hit&autowater.nums==' aw '&stim.num==' stim '&~early']};       % All hits, 2afc, stim on/off, not early
params.modecondition(6) = {['hit&autowater.nums==1&stim.num==' stim '&~early']};        % All hits, aw on, stim on/off, not early

%% SET METADATA FROM ALL RELEVANT SESSIONS/ANIMALS
meta = [];
meta = loadJEB4_ALMVideo(meta);
meta = loadJEB5_ALMVideo(meta);
meta = loadJEB6_ALMVideo(meta);
meta = loadJEB7_ALMVideo(meta);
% meta = loadEKH1_ALMVideo(meta);
meta = loadEKH3_ALMVideo(meta);
meta = loadJGR2_ALMVideo(meta);
meta = loadJGR3_ALMVideo(meta);

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
%%
for gg = 1:length(meta)         % For all loaded sessions...
    ff = figure(gg);
    ff.WindowState = 'maximized';
    obj = objs{gg};
    met = meta(gg);

    anm = obj.pth.anm;                  % Animal name
    date = obj.pth.dt;                  % Session date
    probenum = string(met.probe);       % Which probe was used

    clear rez; clear removeEarly, clear reg

    rez.time = objs{1}.time;
    rez.condition = objs{1}.condition;
    rez.alignEvent = params.alignEvent;

    % Find motion energy for each trial and determine whether or not
    % animal is moving at all points
    edges = 0:0.005:5.5;
    conditions = {1,2};
    [met,mov,me] = assignEarlyTrials(obj,met,params);
    [MEinterp,MEraw] = findInterpME(edges,conditions, met,mov,me);

    % Find all modes
    allModes = calcAllModes(obj,met,rez,params,'no');

    % Project single trials onto choice mode
    cd = allModes.choice_mode;
    conditions = {1,2};
    latent = getTrialLatents(obj,cd,conditions,met);

    % Which conditions to project onto the modes
    conditions = [1,2];         % Left and right 2AFC hits (not early)
    smooth = 61;
    allModes = getChoiceModeProjection(obj,allModes,smooth,conditions);

    % Plot the R and L hits onto the two choice modes (found w/ and w/o
    % early move trials)
    if ~isempty(allModes.latentChoice)
        subplot(3,2,2)
        colors = {[0 0 1],[1 0 0],[0.5 0.5 1],[1 0.5 0.5]};
        lw = 2;
        plot(rez.time,allModes.latentChoice{1},'Color',colors{1},'LineWidth',2)
        hold on
        plot(rez.time,allModes.latentChoice{2},'Color',colors{2},'LineWidth',2)
        legend('Right','Left','Location','best')
        title('Delay CD','FontSize',14)
        xlabel('Time since go-cue (s)','FontSize',13)
        ylabel('Choice Mode (a.u.)','FontSize',13)
        xlim([-2.5 2.5])
    end

    % PANEL A: Plot prob of movement for both trial types
    subplot(3,2,1)
    conditions = {1,2};
    colors = {[0 0 1],[1 0 0],[0.5 0.5 1],[1 0.5 0.5]};
    plotMoveProb_SessAvg(met,conditions,colors,mov,me)
    legend('Right','Left','Location','best')

    % Combining choice mode information across conditions into a single struct
    lat_choice = [];
    movement_int = [];
    movement = [];
    for c = 1:numel(conditions)
        lat_choice = [lat_choice,latent{c}];
        movement_int = [movement_int,MEinterp{c}];
        movement = [movement,MEraw{c}];
    end

    % Find the average jaw velocity during specified time points (on each
    % trial)
    startix = find(taxis>=-0.4, 1, 'first');
    stopix = find(taxis<=-0.05, 1, 'last');
    val = mean(movement_int(startix:stopix, :), 1,'omitnan');
    % val = lastlick;

    % Sort the average jaw velocities in descending order and save the trial
    % order
    [~, ix] = sort(val, 'descend');

    Ntrials = numel(ix);
    Ngroups = 5;                                % Number of groups that you want to partition trials into
    trialsPerGroup = floor(Ntrials/Ngroups);    % Number of trials to include in each group
    num = [1:Ntrials];
    group = ceil(num/trialsPerGroup);            % Assign a group to each sorted trial
    group(group>Ngroups) = Ngroups;             % Any group number that is above the Ngroups, change it to the last group number


    groupsToPlot = [1:5];                       % Which of the groups do you want to plot?
    groupslegend = cell(1,numel(groupsToPlot));

    for i = 1:numel(groupsToPlot)
        temp = groupsToPlot(i);
        if temp == 5
            groupslegend{i} = '5 - Low ME group';
        elseif temp == 1
            groupslegend{i} = '1 - High ME group';
        else
            groupslegend{i} = num2str(temp);
        end
    end

    % Plot a heatmap of all of the all of the single trial choice latents, sorted according to the average jaw velocity
    subplot(3,2,4); imagesc(taxis,1:numel(ix),lat_choice(:, ix)');
    xlabel('Time since go-cue (s)')
    ylabel('Sorted trials')
    title('Single trial choice mode latents')
    colorbar

    % Heatmap of jaw velocity, sorted accordingly
    subplot(3,2,3); imagesc(taxis,1:numel(ix),movement_int(:, ix)');
    xlabel('Time since go-cue (s)')
    ylabel('Sorted trials')
    title('Single trial motion energy')
    colorbar

    % FIGURE 1, Panel A: Plot single trial choice latents
    ax1 = subplot(3,2,5); hold on;
    clr = colormap(ax1,jet(Ngroups));               % Generate Ngroups number of colors from the specified colormap (Ngroups x 3 struct where the 3 is an RGB value)


    for i = 1:Ntrials                 % For all trials...
        trix = ix(i);
        c = clr(group(i), :);                           % Find which group this trial belongs to and assign the color associated with that group
        if ismember(group(i), groupsToPlot)             % If the trial is a member of the groups that you want to plot...
            plot(taxis,medfilt1(lat_choice(:, trix), 151), 'Color', c, 'LineWidth', 1.5);        % Plot the single trial projections onto the choice mode. Colored according to what the jaw velocity was on that trial
        end
    end
    xlabel('Time since go-cue (s)')
    ylabel('Choice mode (a.u.)')
    title('Single trial choice mode latents')
    hold off;

    % FIGURE 1, Panel B: Plot average choice latents for each group
    ax2 = subplot(3,2,6); hold on;
    clr = colormap(ax2,jet(Ngroups));
    for i = 1:Ngroups                 % For all groups of trials...
        c = clr(i, :);                                  % Find the color associated with the current group
        trix = ix(find(group==i));
        if ismember(i, groupsToPlot)                    % If current group is one that you want to plot...
            ts = mean(medfilt1(lat_choice(:, trix), 25), 2,'omitnan');        % Find the average projection onto the choice mode for that group of trials
            plot(taxis,ts, 'Color', c, 'LineWidth', 3);
        end
    end
    legend(groupslegend,'Location','best')
    xlabel('Time since go-cue (s)')
    ylabel('Choice mode (a.u.)')
    title('Avg choice mode latents for each group')
    hold off;

    sesstitle = strcat(anm,date,' ;  ','Probe ',probenum,'ME_LateDelay');  % Name/title for session
    sgtitle(sesstitle,'FontSize',16)

    if strcmp(toSave,'yes')
        saveas(gcf,fullfile(outputdir,sesstitle),'jpeg')
        close all
    end

    %close all;
end
%%
%
%
%
% jv = jaw;
% jv(jv>10) = 10;
% x = 1:size(jaw, 1);
% figure; hold on;
% for i = 1:size(lat_choice, 2)
%    scatter(x, medfilt1(lat_choice(:, i), 75),  50, jv(:, i) , '.');
%
% end
%
%
% lc = medfilt1(lat_choice, 151);
% figure; plot(jv(:), lc(:), '.');
%
% jvm = nanmean(jaw(60:120, :), 1);
% [~, ix] = sort(jvm, 'descend');

% % Get the timing of the last lick in each trial
% lastlick = zeros(size(jaw, 2), 1);      % Trials x 1
% for i = 1:size(jaw, 2)                  % For all trials...
%     licks = sort([obj.bp.ev.lickR{i} obj.bp.ev.lickL{i}]);      % Find all of the licks that occurred on that trial and sort them in time order
%     if isempty(licks)
%         lastlick(i) = 0;
%     else
%         lastlick(i) = licks(end);                                   % Find the time of the last lick in the trial
%     end
% end

