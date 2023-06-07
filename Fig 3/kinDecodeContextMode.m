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
% meta = loadJEB13_ALMVideo(meta,datapth);
% meta = loadJEB14_ALMVideo(meta,datapth); % selectivity in ME % go cue is at 2.3 instead of 2.5 like all other sessions??
% meta = loadJEB15_ALMVideo(meta,datapth);
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

%% CODING DIMENSIONS

clearvars -except obj meta params sel_corr_mat datapth me kin

% % 2afc (early, late, go)
% cond2use = [2 3]; % left hit, right hit
% cond2proj = [2 3 4 5 6 7 8 9];
% rez_2afc = getCodingDimensions_2afc(obj,params,cond2use,cond2proj, 1);

% % aw (context mode)
cond2use = [6 7]; % hit 2afc, hit aw
cond2proj = [2 3 4 5 6 7 8 9];
rez_aw = getCodingDimensions_aw(obj,params,cond2use,cond2proj);

%allrez = concatRezAcrossSessions(rez_2afc,rez_aw);

%% CONTEXT MODE ON SINGLE TRIALS

% cond2use = [6 7]; % hit 2afc, hit aw
for sessix = 1:numel(meta)
    Wcontext = rez_aw(sessix).cd_mode_orth;
    trialdat{sessix} = tensorprod(obj(sessix).trialdat,Wcontext,2,1); % (time,trials)
end

%% DECODING PARAMETERS

% input data = neural data (time*trials,neurons)
% output data = kin data   (time*trials,kin feats)

par.pre=6; % time bins prior to output used for decoding
par.post=0; % time bins after output used for decoding
par.dt = params(1).dt; % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using. Set params.dt and pre/post appropriately for you analysis
par.post_s = par.post .* params(1).dt;

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
par.feats = kin(1).featLeg;
par.feats = {'tongue','motion','nose','jaw'};
temp = cellfun(@(x) patternMatchCellArray(kin(1).featLeg,{x},'all') , par.feats,'UniformOutput',false);
par.feats = cat(1, temp{:});
% par.feats = {'tongue_angle','tongue_length','motion_energy'};

% trials
par.cond2use = [6 7];

par.regularize = 0; % if 0, linear regression. if 1, ridge regression

%% DECODING

close all

for sessix = 1:numel(meta)
    disp([num2str(sessix) ' / ' num2str(numel(meta))])

    % predict 'trialdat' from par.feats (kinematic features of interest)

    % trials
    par.trials.all = cell2mat(params(sessix).trialid(par.cond2use)');

    nTrials = numel(par.trials.all);
    nTrain = floor(nTrials*par.train);
    par.trials.train = par.trials.all; %randsample(par.trials.all,nTrain,false);
    % par.trials.test = par.trials.all(~ismember(par.trials.all,par.trials.train));
    par.trials.test = [];


    % input data
    par.featix = find(ismember(kin(sessix).featLeg,par.feats));

    X.train = kin(sessix).dat(:,par.trials.train,par.featix); % (time,trials,feats)
    X.size = size(X.train);
    X.train = reshape(X.train, size(X.train,1)*size(X.train,2),size(X.train,3));

    X.test = kin(sessix).dat(:,par.trials.test,par.featix); % (time,trials,feats)
    X.test = permute(X.test,[1 3 2]);
    X.test = reshape(X.test, size(X.test,1)*size(X.test,2),size(X.test,3));

    % reshape train and test data to account for prediction bin size
    X.train = reshapePredictors(X.train,par);
    X.test = reshapePredictors(X.test,par);

    % flatten inputs
    % if you're using a model with recurrence, don't flatten
    X.train = reshape(X.train,size(X.train,1),size(X.train,2)*size(X.train,3));
    X.test = reshape(X.test,size(X.test,1),size(X.test,2)*size(X.test,3));

    % output data
    Y.train = trialdat{sessix}(:,par.trials.train); % (time,trials);
    Y.size = size(Y.train);
    Y.train = reshape(Y.train, size(Y.train,1)*size(Y.train,2),size(Y.train,3));

    Y.test = trialdat{sessix}(:,par.trials.test);
    Y.test = reshape(Y.test, size(Y.test,1)*size(Y.test,2),size(Y.test,3));

    % standardize data
    % standardize both train and test sets using train set statistics
    % can also standardize using specific time points (presample for example)
    X.mu = mean(X.train,1,'omitnan');
    X.sigma = std(X.train,[],1,'omitnan');
    X.train = (X.train - X.mu) ./ X.sigma;
    if ~par.test==0
        X.test = (X.test - X.mu) ./ X.sigma;
    end

    Y.mu = mean(Y.train,1,'omitnan');
    Y.sigma = std(Y.train,[],1,'omitnan');
    Y.train = (Y.train - Y.mu);% ./ Y.sigma;
    if ~par.test==0
        Y.test = (Y.test - Y.mu);% ./ Y.sigma;
    end

    % fill missing values in kinematics
    X.train = fillmissing(X.train,'constant',0);
    Y.train = fillmissing(Y.train,'nearest');
    X.test = fillmissing(X.test,'constant',0);
    Y.test = fillmissing(Y.test,'nearest');

    if par.regularize
        mdl = fitrlinear(X.train,Y.train,'Learner','svm','KFold',par.nFolds,'Regularization','ridge');
    else
        mdl = fitrlinear(X.train,Y.train,'Learner','svm','KFold',par.nFolds);
    end
    pred = kfoldPredict(mdl);

    nPredictors = size(X.train,2);                  % Number of kinematic predictors * binWidth
    loadings = NaN(nPredictors,par.nFolds);         % (# predictors x folds for CV)
    for fol = 1:par.nFolds
        loadings(:,fol) = mdl.Trained{fol}.Beta;
    end
    avgloadings(:,sessix) = mean(loadings,2,'omitnan');  % Average the coefficients for each predictor term across folds; save these for each time point

end
%%
binWidth = par.pre-par.post;
featid =[];
for ff = 1:numel(par.feats)
    temp = cell(1,binWidth);
    for bb = 1:binWidth
        temp{bb} = par.feats{ff};
    end
    featid = [featid,temp];
end
%% Normalize all of the beta coefficients to be a fraction of 1
clear temp
featgroups = {'tongue','nos','jaw','motion_energy'};
epochfns = {'delay'};
totalnormfeats = NaN(length(meta),length(featgroups));
for sessix = 1:length(meta)
    totalbeta = sum(abs(avgloadings(:,sessix)));
    allLoadings(sessix).totalbeta = totalbeta;
    for group = 1:length(featgroups)
        temp = zeros(1,length(featid));
        currgroup = featgroups{group};
        for feat = 1:length(featid)
            currfeat = featid{feat};
            temp(feat) = contains(currfeat,currgroup);
        end
        groupixs = find(temp);
        allLoadings(sessix).(currgroup) = abs(avgloadings(groupixs, sessix));
    end
end

totalnormfeats = NaN(length(meta),length(featgroups));
for group = 1:length(featgroups)
    currgroup = featgroups{group};
    for sessix = 1:length(meta)
        groupLoad = allLoadings(sessix).(currgroup);
        grouptotal = sum(groupLoad,'omitnan');
        grouprelative = grouptotal / allLoadings(sessix).totalbeta;
        grouprelative = grouprelative/(length(groupLoad));
        totalnormfeats(sessix,group) = grouprelative;
    end
end

for sessix = 1:length(meta)
    currsess = totalnormfeats(sessix,:);
    tot = sum(currsess,'omitnan');
    totalnormfeats(sessix,:) = currsess./tot;
end

% Get all ME vals
MEix = find(strcmp(featgroups,'motion_energy'));
MEvals = totalnormfeats(:,MEix);
[~,sortix] = sort(MEvals,'descend');
totalnormfeats = totalnormfeats(sortix,:);
%% Make a stacked bar plot to show across sessions that it varies which feature group contributes most to predicting CDTrialType
bar(totalnormfeats,'stacked')
legend(featgroups,'Location','best')
xlabel('Session #')
ylabel('Sum of beta coefficients for each feature group')

%%
feat2use = {'jaw_yvel_view1','nose_yvel_view1', 'tongue_length'};
featix = NaN(1,length(feat2use));
for f = 1:length(feat2use)
    currfeat = feat2use{f};
    currix = find(strcmp(kin(1).featLeg,currfeat));
    featix(f) = currix;
end

cond2use = [6 7];
condfns = {'2AFC hit','AW hit'};
trix2use = 20;

times.start = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent));
times.startix = find(obj(1).time>times.start,1,'first');
times.stop = mode(obj(1).bp.ev.sample)-mode(obj(1).bp.ev.(params(1).alignEvent))-0.05;
times.stopix = find(obj(1).time<times.stop,1,'last');

del = median(obj(1).bp.ev.delay)-median(obj(1).bp.ev.(params(1).alignEvent));
delix = find(obj(1).time>del,1,'first');
go = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent));
goix = find(obj(1).time<go,1,'last');
resp = median(obj(1).bp.ev.goCue)-median(obj(1).bp.ev.(params(1).alignEvent))+2.5;
respix = find(obj(1).time<resp,1,'last');

sm = 31;
ptiles = [92 95 20];

sortedSess2use = 11;
sessix = sortix(sortedSess2use);


for c = 1:length(cond2use)
    allkin = [];
    figure();
    cond = cond2use(c);
    condtrix = params(sessix).trialid{cond};
    ntrials = length(condtrix);
    randtrix = randsample(condtrix,trix2use);
    for f = 1:length(featix)
        currfeat = featix(f);
        
        % ^ Don't actually want to max normalize because then will be
        % normalizing by an outlier value probably

        % Want to normalize to the 90-99th percentile of values to account
        % for more of the data
   
        if strcmp(feat2use{f},'tongue_length')
            currkin = kin(sessix).dat(times.startix:goix,randtrix,currfeat);
            abskin = abs(currkin);
            normkin = abskin./prctile(abskin(:), ptiles(f));
            normkin(normkin>1) = 1;
        else
            currkin = mySmooth(kin(sessix).dat_std(times.startix:goix,randtrix,currfeat),sm);
            abskin = abs(currkin);
            normkin = abskin./prctile(abskin(:), ptiles(f));
            normkin(normkin>1) = 1;
        end                                                                  % Will end up with values greater than 1 in this case--set these to 1

        allkin = cat(3,allkin,normkin);                                      % Concatenate across features (trials x time x feat)
    end

    allkin = permute(allkin,[2 1 3]);                                        % (time x trials x feat/RGB)
    RI = imref2d(size(allkin));
    RI.XWorldLimits = [0 3];
    RI.YWorldLimits = [2 5];
    IMref = imshow(allkin, RI,'InitialMagnification','fit');
    title(['RGB = ' feat2use '; Sorted session ' num2str(sortedSess2use) ' ; ' condfns{c}])
end
%%
sortedSess2use = 6;
sess2use = sortix(sortedSess2use);

% Plot all features for the same trial in one subplot

colors = {[1 0 0],[0 1 0],[0 0 1]};
nTrixPlot = 9;
offset = 3;
sm = 20;
for sessix = sess2use
    for c = 1:length(cond2use)
        figure();
        condtrix = params(sessix).trialid{cond2use(c)};
        trix = randsample(condtrix,nTrixPlot);
        for tt = 1:nTrixPlot
            subplot(3,3,tt)
            for feat = 1:length(featix)
                if strcmp(feat2use{feat},'tongue_length')
                    condkin = kin(sessix).dat(:,trix(tt),featix(feat));
                    condkin = condkin./prctile(condkin(:), 1);
%                     condkin(condkin>1) = 1;
                    toplot = offset*feat+abs(condkin(:,:,:));
                else
                    condkin = kin(sessix).dat_std(:,trix(tt),:);
                    toplot = offset*feat+abs(mySmooth(condkin(:,:,featix(feat)),sm));
                end
                plot(obj(sessix).time,toplot,'LineWidth',2,'Color',colors{feat}); hold on;
            end
            title(feat2use(feat))
            xlabel('Time from go cue (s)')
            ylim([2 14])
            xlim([-2.5 0])
            xline(-0.9,'k--','LineWidth',1)
            xline(-2.2,'k--','LineWidth',1)
        end
        sgtitle(['Sorted session ' num2str(sortedSess2use) ' ; ' condfns{c}])
    end
end