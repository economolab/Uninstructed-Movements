% DECODING CDlate FROM ALL KINEMATIC FEATURES
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
meta = loadJEB13_ALMVideo(meta,datapth);
meta = loadJEB6_ALMVideo(meta,datapth);
meta = loadJEB7_ALMVideo(meta,datapth);
meta = loadEKH1_ALMVideo(meta,datapth);
meta = loadEKH3_ALMVideo(meta,datapth);
meta = loadJGR2_ALMVideo(meta,datapth);
meta = loadJGR3_ALMVideo(meta,datapth);
meta = loadJEB14_ALMVideo(meta,datapth);
meta = loadJEB15_ALMVideo(meta,datapth);
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
    disp(['Loading ME for session ' num2str(sessix)])
    me(sessix) = loadMotionEnergy(obj(sessix), meta(sessix), params(sessix), datapth);
end
%% Calculate all CDs and find single trial projections
clearvars -except obj meta params me sav kin

disp('----Calculating coding dimensions----')
cond2use = [2 3 6 7]; % right hits, left hits (corresponding to PARAMS.CONDITION)
inclramp = 'yes';
rampcond = 8;
cond2proj = 2:7;  % right hits, left hits, right miss, left miss, right no, left no (corresponding to null/potent psths in rez)
cond2use_trialdat = [2 3]; % for calculating selectivity explained in full neural pop
regr = getCodingDimensions_2afc(obj,params,cond2use,rampcond,cond2proj, inclramp);

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
%% Predict CDTrialType from DLC features
clearvars -except datapth kin me meta obj params regr nSessions exsess
clc;

% input data = neural data (time*trials,neurons)
% output data = kin data   (time*trials,kin feats)
par.pre=5; % time bins prior to output used for decoding
par.post=0; % time bins after output used for decoding
par.dt = params(1).dt; % moving time bin
par.pre_s = par.pre .* params(1).dt; % dt, pre_s, and post_s just here to know how much time you're using. Set params.dt and pre/post appropriately for you analysis
par.post_s = par.post .* params(1).dt;


trialstart = mode(obj(1).bp.ev.bitStart)-mode(obj(1).bp.ev.(params(1).alignEvent)+0.4);
start = find(obj(1).time>trialstart,1,'first');
go = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent));
stop = find(obj(1).time>go,1,'first');

par.timerange = start:stop;

% cross val folds
par.nFolds = 4;

% data sets
par.train = 1; % fraction of trials (just using cross-val here, matlab's kFoldPredict uses held out data for testing)
par.test = 1 - par.train;

% feature to use to decode
% par.feats = kin(1).featLeg;
par.feats = {'motion','nos','jaw'};
temp = cellfun(@(x) patternMatchCellArray(kin(1).featLeg,{x},'all') , par.feats,'UniformOutput',false);
par.feats = cat(1, temp{:});
% par.feats = {'tongue_angle','tongue_length','motion_energy'};

% trials
par.cond2use = [2 3];

par.regularize = 1; % if 0, linear regression. if 1, ridge regression
regtype = 'Ridge';
%% DECODING
close all

[avgloadings, trueVals, modelpred] = doDLC_CDChoiceDecoding(meta, par, kin, regr, params);

disp('---FINISHED DECODING FOR ALL SESSIONS---')
clearvars -except datapth kin me meta obj params regr nSessions exsess modelpred trueVals par avgloadings
%% Show how the movements of different body parts are driving the potent space
binWidth = par.pre-par.post;
featid =[];
for ff = 1:numel(par.feats)
    temp = cell(1,binWidth);
    for bb = 1:binWidth
        temp{bb} = par.feats{ff};
    end
    featid = [featid,temp];
end

del = mode(obj(1).bp.ev.delay)-mode(obj(1).bp.ev.(params(1).alignEvent)+0.4);
start = find(obj(1).time>del,1,'first');
resp = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent))-0.05;
stop = find(obj(1).time<resp,1,'last');
resp2 = mode(obj(1).bp.ev.goCue)-mode(obj(1).bp.ev.(params(1).alignEvent))+2;
stop2 = find(obj(1).time<resp2,1,'last');
%% Normalize all of the beta coefficients to be a fraction of 1
clear temp
featgroups = {'nos','jaw','motion_energy'};
epochfns = {'delay'};

allLoadings = getAllLoadingsForFeat(meta, avgloadings,featgroups, featid);

totalnormfeats = normalizeFeatureLoadings(meta, featgroups, allLoadings);

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

clearvars -except kin regr meta obj params modelpred trueVals totalnormfeats sortix avgloadings par
%%
feat2use = {'jaw_yvel_view1','nose_yvel_view1', 'motion_energy'};
featix = NaN(1,length(feat2use));
for f = 1:length(feat2use)
    currfeat = feat2use{f};
    currix = find(strcmp(kin(1).featLeg,currfeat));
    featix(f) = currix;
end

cond2use = [2 3];
condfns = {'R hit','Lhit'};
trix2use = 40;

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
ptiles = [94 98 94];

sortedSess2use = [1 25];

for sessix = 9 %ss = 1:length(sortedSess2use)
    %sessix = sortix(sortedSess2use(ss));

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

            if strcmp(feat2use{f},'motion_energy')
                currkin = mySmooth(kin(sessix).dat(times.startix:goix,randtrix,currfeat),sm);
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
%         MEix = find(strcmp(feat2use,'motion_energy'));
%         currME = squeeze(allkin(:,:,MEix));
%         delME = currME(:,(size(currME,2)/2):end);
%         delME = mean(delME,2,'omitnan');
%         [~,sortix] = sort(delME,'descend');
%         allkin(:,:,MEix) = currME(sortix,:);

        RI = imref2d(size(allkin));
        RI.XWorldLimits = [0 3];
        RI.YWorldLimits = [2 5];
        IMref = imshow(allkin, RI,'InitialMagnification','fit');
%         title(['RGB = ' feat2use '; Sorted session ' num2str(sortedSess2use(ss)) ' ; ' condfns{c}])
        title(['RGB = ' feat2use '; Sorted session ' num2str(sessix) ' ; ' condfns{c}])
    end
end
%%
sess2use = sessix;

% Plot all features for the same trial in one subplot

colors = {[1 0 0],[0 1 0],[0 0 1]};
nTrixPlot = 9;
offset = 3;
sm = 31;
for sessix = sess2use
    for c = 1:length(cond2use)
        figure();
        condtrix = params(sessix).trialid{cond2use(c)};
        trix = randsample(condtrix,nTrixPlot);
        for tt = 1:nTrixPlot
            subplot(3,3,tt)
            for feat = 1:length(featix)
                if strcmp(feat2use{feat},'motion_energy')
                    condkin = kin(sessix).dat(:,trix(tt),featix(feat));
                    condkin = condkin./prctile(condkin(:), 60);
%                     condkin(condkin>1) = 1;
                    toplot = offset*feat+abs(mySmooth(condkin(:,:,:),5));
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