% DECODING CDlate FROM ALL KINEMATIC FEATURES
clear,clc,close all

% add paths for data loading scripts, all fig funcs, and utils
utilspth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code\Uninstructed-Movements\NullPotent_Context';
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
% rmpath(genpath(fullfile(utilspth,'fig1/')));
% rmpath(genpath(fullfile(utilspth,'mc_stim/')));

% add paths for figure specific functions
addpath(genpath(pwd))

%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                             % all trials
params.condition(end+1) = {'R&hit&~stim.enable&~autowater'};             % R 2AFC hits, no stim
params.condition(end+1) = {'L&hit&~stim.enable&~autowater'};             % L 2AFC hits, no stim
params.condition(end+1) = {'R&miss&~stim.enable&~autowater'};            % R error 2AFC, no stim, aw off
params.condition(end+1) = {'L&miss&~stim.enable&~autowater'};            % L error 2AFC, no stim

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance

params.N_varToExplain = 80; % keep num dims that explains this much variance in neural data (when doing n/p)

params.advance_movement = 0;
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
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);

% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

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
%% Calculate all CDs and find single trial projections
clearvars -except obj meta params me sav kin

disp('----Calculating coding dimensions----')
cond2use = [2,3];
cond2proj = [2,3];
regr = getCodingDimensions(obj,params,cond2use,cond2proj);

disp('----Projecting single trials onto CDlate----')
cd = 'late';
regr = getSingleTrialProjs(regr,obj,cd);
%% Load kinematic data
nSessions = numel(meta);
for sessix = 1:numel(meta)
    message = strcat('----Getting kinematic data for session',{' '},num2str(sessix), {' '},'out of',{' '},num2str(nSessions),'----');
    disp(message)
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% CDLate from DLC features

clearvars -except datapth kin me meta obj params regr nSessions

% params
rez.nFolds = 4;                                     % number of iterations (bootstrap)

rez.binSize = 30;                                   % Bin size to decode over (in milliseconds)
difTime = params(1).dt*1000;                        % Convert time-difference from s to ms
rez.dt = floor(rez.binSize / difTime);              % How many samples confer your bin size
rez.tm = obj(1).time(1:rez.dt:numel(obj(1).time));  % Create a new time axis over which to decode (based on the bin size that you want)
rez.numT = numel(rez.tm);                           % Number of time-points

rez.train = 0.5;                                      % fraction of trials to use for training (1-train for testing)

% match number of right and left hits, and right and left misses
cond2use = 2:5;
hitcond = [1 2];
misscond = [3 4];

% A lot of features have redundant measurements--group all measurements
% for the same feature together
featGroups = {{'tongue'},...
    {'jaw','trident'},...
    {'nose','nostril'},...
    {'paw'},...
    {'motion_energy'}};


% True Values of CDlate for each session
trueVals.Rhit = cell(nSessions,1);
trueVals.Lhit = cell(nSessions,1);

% Model prediction of CDlate for each session, predicted by each feature
modelpred.Rhit = cell(nSessions,numel(featGroups));
modelpred.Lhit = cell(nSessions,numel(featGroups));
for ifeat = 1:numel(featGroups)
    % Specify which body part you want to decode with
    disp(['----Feature Group ' num2str(ifeat) '/' num2str(numel(featGroups)) '----'])
    %     rez.feats2use = kin(1).featLeg;
    % rez.feats2use = {'jaw_ydisp_view1'};
    % rez.feats2use = {'motion_energy'};
    % rez.feats2use = {'view2'};
    rez.feats2use = featGroups{ifeat};
    if size(rez.feats2use,1) == 1
        rez.feats2use = rez.feats2use';
    end


    % Find which indices in the feature legend correspond to the specified body part 
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
    rez.featix = find(mask);                    % Indices in the feature legend that correspond to specified body part


    % Do the decoding for each session
    for sessix = 1:numel(obj)
        disp(['Decoding session ' num2str(sessix) ' / ' num2str(numel(obj))])

        % Getting the proper number of trials
        trials_cond = params(sessix).trialid(cond2use);                                         % Correct trials for each condition you use   

        minHitTrials = cellfun(@(x) numel(x),trials_cond(hitcond), 'UniformOutput',false);      % Find the number of trials for all hit conditions (i.e. R and L) -- apply the function (numel) to each element in trialcond    
        nhits = min(cell2mat(minHitTrials));                                                    % Which condition had the fewer number of trials 

        minMissTrials = cellfun(@(x) numel(x),trials_cond(misscond), 'UniformOutput',false);    % Find which condition had fewer number of miss trials
        nmiss = min(cell2mat(minMissTrials));


        trials_hit = cellfun(@(x) randsample(x,nhits), trials_cond(hitcond), 'UniformOutput', false);       % From each of the condition trials, sample the proper number of trials for the training set
        trialsHit = cell2mat(trials_hit);                                                                   % Convert cell array into a double
        trialsHit = trialsHit(:);                                                                           % Stack the columns of trials on top of each other

        trials_miss = cellfun(@(x) randsample(x,nmiss), trials_cond(misscond), 'UniformOutput', false);
        trialsMiss = cell2mat(trials_miss);
        trialsMiss = trialsMiss(:);

        trials.all = [trialsHit ; trialsMiss];                                                              % Take all hit and miss trials that you are using for training and put them in one column

        % OUTPUT: What you are trying to predict--CDlate for all of the train trials
        % Y = time x trials
        Y = [];
        for t = 1:length(trials.all)
            currtrix = trials.all(t);
            Y = [Y,regr(sessix).singleProj(:,currtrix)];
        end

        % INPUT: What you are using to predict--kinematic measures for the given body part for all of the trian trials
        % X = time x trials x num measures for body part
        X = kin(sessix).dat(:,trials.all,rez.featix);
        % fill missing values
        for featix = 1:size(X,3)
            X(:,:,featix) = fillmissing(X(:,:,featix),"constant",0);
        end

        % Train/Test Split

        % Randomly sample (without replacement) the specified number of
        % trials from the data--based on the fraction you specified in
        % 'rez.train'
        [trials.train,trials.trainidx] = datasample(trials.all,round(numel(trials.all)*rez.train),'Replace',false);     % Trial numbers for training set, indices within trials.all that correspond to the Training trials
        trials.testidx = find(~ismember(trials.all,trials.train));                                                      % Find the indices that are not in the training set (i.e. the test set)
        trials.test = trials.all(trials.testidx);                                                                       % Set those trials as the test set

        % Get the predictor and regressor data for the train/test split 
        in.train.y = Y(:,trials.trainidx);
        in.test.y  = Y(:,trials.testidx);
        in.train.X = X(:,trials.trainidx,:);
        in.test.X  = X(:,trials.testidx,:);

        % decoding
        pred = DLC_CD_Decoder(in,rez);
        pred = pred';
        
        % Which of the test trials used were R and L hits?
        trials.RHit.TestIX = ismember(trials.test,trials_hit{1});       % Logical array for all of the test trials indicating whether they were a right hit or not
        trials.RHit.Test = trials.test(trials.RHit.TestIX);             % Get the trial numbers that are a R hit and were used in the test
        trials.LHit.TestIX = ismember(trials.test,trials_hit{2});       % Logical array for all of the test trials indicating whether they were a left hit or not
        trials.LHit.Test = trials.test(trials.LHit.TestIX);             % Get the trial numbers that are a L hit and were used in the test
        
        % Divide test data used and model prediction into R and L
        if ifeat==1
            trueVals.Rhit{sessix} = in.test.y(:,trials.RHit.TestIX);
            trueVals.Lhit{sessix} = in.test.y(:,trials.LHit.TestIX);
        end
        modelpred.Rhit{sessix,ifeat} = pred(:,trials.RHit.TestIX);
        modelpred.Lhit{sessix,ifeat} = pred(:,trials.LHit.TestIX);
       
        % %     % shuffle labels for a 'null' distribution
        % %
        % %
        % %     Y = randsample(Y,numel(Y));
        % %
        % %     % train/test split
        % %
        % %     in.train.y = Y(trials.trainidx);
        % %     in.test.y  = Y(trials.testidx);
        % %
        % %     acc_shuffled(:,sessix) = DLC_ChoiceDecoder(in,rez,trials);
    end

end
%% Plot an example session of CDlate prediction vs true value
sessix = 1;
featGroup = 2;
avgCD.Rhit.true = mean(trueVals.Rhit{sessix},2);
avgCD.Lhit.true = mean(trueVals.Lhit{sessix},2);

avgCD.Rhit.pred = mean(modelpred.Rhit{sessix,featGroup},2);
avgCD.Lhit.pred = mean(modelpred.Lhit{sessix,featGroup},2);


figure()
subplot(2,1,1)
plot(obj(1).time,trueVals.Rhit{sessix},'Color','blue','LineWidth',2); hold on;
plot(rez.tm(1:end-1),modelpred.Rhit{sessix,featGroup},'Color',[0.5 0.5 1],'LineWidth',2);

subplot(2,1,2)
plot(obj(1).time,trueVals.Lhit{sessix},'Color','red','LineWidth',2); hold on;
plot(rez.tm(1:end-1),modelpred.Lhit{sessix,featGroup},'Color',[1 0.5 0.5],'LineWidth',2);
ylabel('a.u.')
xlabel('Time since go-cue (s)')

figure()
plot(obj(1).time,avgCD.Rhit.true,'Color','blue','LineWidth',2); hold on;
plot(rez.tm(1:end-1),avgCD.Rhit.pred,'Color',[0.5 0.5 1],'LineStyle','--','LineWidth',2);
plot(obj(1).time,avgCD.Lhit.true,'Color','red','LineWidth',2); hold on;
plot(rez.tm(1:end-1),avgCD.Lhit.pred,'Color',[1 0.5 0.5],'LineStyle','--','LineWidth',2);
legend('R hit true','R hit predicted','L hit true','L hit predicted')
ylabel('a.u.')
xlabel('Time since go-cue (s)')