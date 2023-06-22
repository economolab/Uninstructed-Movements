clear,close all

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
addpath(genpath(fullfile(figpth,'Hazarded Delay')));
figpth = [basepth  '\Uninstructed-Movements\functions'];
addpath(genpath(fullfile(figpth,'HazDel funcs')));
%% PARAMETERS
params.alignEvent          = 'delay'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early'};             % right hits, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};            % error right, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};             % left hits, no stim, aw off
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};            % error left, no stim, aw off

params.tmin = -2.4;
params.tmax = 2.5;
params.dt = 1/50;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


% params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
%     {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance


params.advance_movement = 0.0;

% Haz delay params
params.delay(1) = 0.3000;
params.delay(2) = 0.6000;
params.delay(3) = 1.2000;
params.delay(4) = 1.8000;
params.delay(5) = 2.4000;
params.delay(6) = 3.6000;
%% SPECIFY DATA TO LOAD

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end
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
%%
% ------------------------------------------
% -- Motion Energy --
% me (struct array) - one entry per session
% ------------------------------------------
for sessix = 1:numel(meta)
    disp(['Loading ME for session ' num2str(sessix)])
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% For each session--Trials separated by delay length
% Get the trialIDs corresponding to each delay length
% Find the PSTH for R and L trials of each delay length
% Find the jaw velocities for R and L trials of each delay length
conditions = [2,3];
condfns = {'Rhit','Lhit'};
for sessix = 1:length(meta)
    del(sessix).delaylen = obj(sessix).bp.ev.goCue - obj(sessix).bp.ev.delay;       % Find the delay length for all trials
    del(sessix).del_trialid = getDelayTrix(params(sessix),conditions,del(sessix));  % Group the trials in each condition based on their delay length
    del(sessix).delPSTH = getPSTHbyDel(params(sessix),del(sessix),obj(sessix), condfns, conditions);             % Get avg PSTH for each delay length
end
%% choice decoding from dlc features

clearvars -except datapth kin me meta obj params

% params
rez.nFolds = 4; % number of iterations (bootstrap)

rez.binSize = 75; % ms
rez.dt = floor(rez.binSize / (params(1).dt*1000)); % samples
rez.tm = obj(1).time(1:rez.dt:numel(obj(1).time));
rez.numT = numel(rez.tm);

rez.train = 1; % fraction of trials to use for training (1-train for testing)

rez.nShuffles = 2;

% match number of right and left hits, and right and left misses
cond2use = 2:5;
hitcond = [1 3];
misscond = [2 4];

% featGroups = {{'tongue'},...
%     {'jaw','trident'},...
%     {'nose','nostril'},...
%     {'paw'},...
%     {'motion_energy'}};

featGroups = {'all'};

for ifeat = 1:numel(featGroups)
    disp(['Feature Group ' num2str(ifeat) '/' num2str(numel(featGroups))])
    %     rez.feats2use = kin(1).featLeg;
    % rez.feats2use = {'jaw_ydisp_view1'};
    % rez.feats2use = {'motion_energy'};
    % rez.feats2use = {'view2'};


    rez.feats2use = featGroups{ifeat};

    if strcmpi(rez.feats2use,'all')
        rez.featix = 1:numel(kin(1).featLeg);
    else

        if size(rez.feats2use,1) == 1
            rez.feats2use = rez.feats2use';
        end

        [~,mask] = patternMatchCellArray(kin(1).featLeg,rez.feats2use,'any');
        if sum(mask) == 0
            [~,mask] = patternMatchCellArray(kin(1).featLeg,rez.feats2use,'all');
        end
        if sum(mask) == 0
            use = zeros(size(kin(1).featLeg));
            for i = 1:numel(kin(1).featLeg)
                for j = 1:numel(rez.feats2use)
                    if contains(kin(1).featLeg{i},rez.feats2use{j})
                        use(i) = 1;
                    end
                end
            end
            mask = logical(use);
        end
        if sum(mask) == 0
            error('didnt find features to use')
        end
        rez.featix = find(mask);

    end



    for sessix = 1:numel(obj)
        disp(['Decoding session ' num2str(sessix) ' / ' num2str(numel(obj))])

        % trials

        trials_cond = params(sessix).trialid(cond2use);

        minHitTrials = cellfun(@(x) numel(x),trials_cond(hitcond), 'UniformOutput',false);
        nhits = min(cell2mat(minHitTrials));

        minMissTrials = cellfun(@(x) numel(x),trials_cond(misscond), 'UniformOutput',false);
        nmiss = min(cell2mat(minMissTrials));


        trials_hit = cellfun(@(x) randsample(x,nhits), trials_cond(hitcond), 'UniformOutput', false);
        trialsHit = cell2mat(trials_hit);
        trialsHit = trialsHit(:);

        trials_miss = cellfun(@(x) randsample(x,nmiss), trials_cond(misscond), 'UniformOutput', false);
        trialsMiss = cell2mat(trials_miss);
        trialsMiss = trialsMiss(:);

%         trials.all = [trialsHit ; trialsMiss];
        trials.all = trialsHit;

        % labels (1 for right choice, 0 for left choice)
        Y = [ones(nhits,1) ; -ones(nhits,1) ; ones(nmiss,1) ; -ones(nmiss,1)]; % right hits, left hits, left miss, right miss

        % input
        % use all features
        X = kin(sessix).dat(:,trials.all,rez.featix);

        % fill missing values
        for featix = 1:size(X,3)
            %             X(:,:,featix) = fillmissing(X(:,:,featix),"constant",0);
        end

        % train/test split

        [trials.train,trials.trainidx] = datasample(trials.all,round(numel(trials.all)*rez.train),'Replace',false);
        trials.testidx = find(~ismember(trials.all,trials.train));
        trials.test = trials.all(trials.testidx);

        in.train.y = Y(trials.trainidx);
        in.test.y  = Y(trials.testidx);
        in.train.X = X(:,trials.trainidx,:);
        in.test.X  = X(:,trials.testidx,:);

        % decoding

        acc(:,sessix,ifeat) = DLC_ChoiceDecoder(in,rez,trials); % (time,sessions,features)


        % shuffle labels for a 'null' distribution


        Y = randsample(Y,numel(Y));

        % train/test split

        in.train.y = Y(trials.trainidx);
        in.test.y  = Y(trials.testidx);
        
        for ishuf = 1:rez.nShuffles
            acc_shuf(:,sessix,ishuf,ifeat) = DLC_ChoiceDecoder(in,rez,trials); % (times,session,ishuf)
        end


    end

end

acc_shuf_ = reshape(acc_shuf,size(acc_shuf,1),size(acc_shuf,2)*size(acc_shuf,3)); % (time,session*ishuf)






