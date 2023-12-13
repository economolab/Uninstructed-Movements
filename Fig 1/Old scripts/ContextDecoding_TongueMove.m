clear,clc,close all

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
figpth = [basepth  '\Uninstructed-Movements\Fig 1'];
addpath(genpath(fullfile(figpth,'funcs')));
%% PARAMETERS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'(hit|miss|no)'};                                              % all trials       (1)
params.condition(end+1) = {'R&hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};      % right hits, 2afc (2)
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};      % left hit, 2afc   (3)
params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};     % right miss, 2afc (4)
params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};     % left miss, 2afc  (5)
params.condition(end+1) = {'hit&~stim.enable&~autowater&~early&((1:Ntrials)>20)'};        % 2afc hits        (6)
params.condition(end+1) = {'hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};         % aw hits          (7)
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % right hits, aw   (8)
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early&((1:Ntrials)>20)'};       % left hits, aw    (9)


params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/100;

% smooth with causal gaussian kernel
params.smooth = 15;
params.bctype = 'reflect'; % reflect, zeropad, none

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality


params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw','top_paw','bottom_paw','top_nostril','bottom_nostril'}};

params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this much variance


params.advance_movement = 0.0;


%% SPECIFY DATA TO LOAD

if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

meta = [];

% --- ALM ---
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth); % selectivity in ME
meta = loadEKH1_ALMVideo(meta,datapth); % selectivity in ME
meta = loadEKH3_ALMVideo(meta,datapth); % selectivity in ME
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB19_ALMVideo(meta,datapth);


% --- M1TJ ---
% meta = loadJEB14_M1TJVideo(meta,datapth);

params.probe = {meta.probe}; % put probe numbers into params, one entry for element in meta, just so i don't have to change code i've already written


%% LOAD DATA

[obj,params] = loadSessionData(meta,params);

for sessix = 1:numel(meta)
    disp(['Loading kinematics - Session ' num2str(sessix) '/' num2str(numel(meta))])
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
    kin(sessix) = getKinematics(obj(sessix), me(sessix), params(sessix));
end
%% Get tongue length and angle for individual licks for all animals and sessions
% Stored in variable 'tonguefeat': (1 x nSessions)
% Each session will have fields 'tongue_length' and 'tongue_angle'
% Each field will have subfields 'RAFC', 'LAFC','RAW','LAW'

nLicks = 5;
lickDur = 10;
cond2use = [2,3,8,9];
condfns = {'RAFC','LAFC','RAW','LAW'};
feat2use = {'tongue_angle','tongue_length'};

featix = getFeatIx(feat2use,kin);

% Get times for which you want to look at the tongue
go = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent));
times.goix = find(obj(1).time<go,1,'last');
resp = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent))+2.5;
times.respix = find(obj(1).time<resp,1,'last');

% tonguefeat = extractAllLicks(meta,feat2use,cond2use,params,times,nLicks,lickDur,kin,featix,condfns);
rawtonguekin = getRawTongueKin(meta,feat2use,cond2use,params,times,nLicks,lickDur,kin,featix,condfns);

clearvars -except tonguefeat obj params meta kin me cond2use condfns feat2use nLicks rawtonguekin
%% DECODING PARAMETERS

% input data = tongue kinematic data
% output data = classification of 2AFC or AW

% par.pre=6; % time bins prior to output used for decoding
% par.post=0; % time bins after output used for decoding
% par.dt = params(1).dt; % moving time bin
% par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using. Set params.dt and pre/post appropriately for you analysis
% par.post_s = par.post .* params(1).dt;

% these parameters above are important for decoding accuracy. for actual
% analyses (like a final analysis to be included in a publication),
% you should vary these parameters only if you have a validation
% set that you won't test on until these parameters are set. Otherwise,
% there's a high risk of overfitting

% cross val folds
par.nFolds = 4;

% data sets
par.train = 1; % fraction of trials (just using cross-val here, matlab's kFoldPredict uses held out data for testing)
par.test = 1 - par.train;

% feature to use to decode
par.feats = {'tongue_angle','tongue_length'};
par.licks2use = 1;
par.maxDur = 10;

% iterations
par.nIterations = 100;

%% DECODING
%%% Use fitcvsm with predictor data being the first lick tongue length for
%%% each trial
close all

accuracy = NaN(length(meta),par.nIterations);
shuffaccuracy = NaN(length(meta),par.nIterations);

for sessix = 1:numel(meta)
    disp(['Decoding session ' ' ' num2str(sessix) ' / ' num2str(numel(meta))])

    for ii = 1:par.nIterations
        % predict context from the tongue angle and tongue length of each first lick in a trial

        % trials
        % get the minimum number of trials across all conditions (RAFC, LAFC, RAW, LAW)
        % subsample trials from each condition to make sure all conditions are trained on the same number of trials
        par = getTrialsforDecoding(rawtonguekin(sessix),condfns,par);

        % input data
        X.train = formatTongueRegressors(par,condfns,rawtonguekin(sessix));         % (nTrainTrials * conds,lick length * num of licks * num regressors)

        % output data
        Y.train = formatOutputData(par,condfns);                            % (nTrainTrials * conds, 1)
        shuffTrix = randsample(1:length(Y.train),length(Y.train),'false');
        Y.shuff = Y.train(shuffTrix);

        % fill missing values in kinematics
        X.train = fillmissing(X.train,'constant',0);

        % fit the SVM
        mdl = fitcsvm(X.train,Y.train,'KFold',par.nFolds);
        pred = kfoldPredict(mdl);

        % fit the SVM on data with shuffled labels
        shuffmdl = fitcsvm(X.train,Y.shuff,'KFold',par.nFolds);
        shuffpred = kfoldPredict(shuffmdl);

        accuracy(sessix,ii) = sum(pred==Y.train)/length(Y.train);
        shuffaccuracy(sessix,ii) = sum(shuffpred==Y.shuff)/length(Y.shuff);
    end
end
clearvars -except accuracy shuffaccuracy obj meta params par kin 
%%
avgAcc_persesh = mean(accuracy,2,'omitnan');
shuffAcc_persesh = mean(shuffaccuracy,2,'omitnan');
X = [1,2];
Y = [mean(avgAcc_persesh),mean(shuffAcc_persesh)];
bar(X,Y); hold on
scatter(1,avgAcc_persesh,'filled','black')
scatter(2,shuffAcc_persesh,'filled','black')

xticks([1,2])
xticklabels({'Data', 'Shuffled'})
ylabel('CV Accuracy of SVM decoder')
title('Context decoding from first lick')