clear,clc,close all

addpath(genpath(pwd))

%% SET RUN PARAMS
params.alignEvent          = 'goCue'; % 'jawOnset' 'goCue'  'moveOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)         = {'(hit|miss|no)'};         % all trials
params.condition(end+1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1)     = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
params.condition(end+1)     = {'R&miss&~stim.enable&~autowater&~early'};        % error right, no stim, aw off
params.condition(end+1)     = {'L&miss&~stim.enable&~autowater&~early'};        % error left, no stim, aw off

% specify probe number and areas to load and process data
params.probe(1) = 1;
params.probeArea{1} = 'ALM';

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

%% SET METADATA
% experiment meta data
% meta.datapth = '/Users/Munib/Documents/Economo-Lab/data';
% meta.datapth = '/Volumes/MUNIB_SSD/Experiments';
meta.datapth = 'M:\Economo-Lab\data';
meta.anm = 'JGR2'; % 'JEB7'  'EKH3'
meta.date = '2021-11-17'; % '2021-04-29'   '2021-08-11
meta.datafn = findDataFn(meta);

%% LOAD AND PROCESS DATA
dat = load(fullfile(meta.datapth,'DataObjects',meta.anm,meta.datafn));
obj = dat.obj;

for prbix = 1:numel(params.probe)
    prbnum = params.probe(prbix);
    disp('______________________________________________________')
    disp(['Processing data for probe ' num2str(prbix)])
    disp(' ')
    [params,obj] = processData(obj,params,prbnum);
end

% if only one probe, clean up so all previous code works
if numel(params.probe)==1
    [obj,params] = oneProbeTrim(obj,params);
end

disp(' ')
disp('DATA LOADED AND PROCESSED')
disp(' ')

%% CURATE DATA FOR GPFA

trials = 1:obj.bp.Ntrials; % all trials
% trials = sort([params.trialid{1} ; params.trialid{2}]); % right and left hits

trialspikes = obj.trialspikes(:,:,trials);

clear dat

for i = numel(trials):-1:1
    dat(i).trialId = trials(i);
    dat(i).spikes  = trialspikes(:,:,i)'; % (neurons,time)
end

%% GPFA FITTING

gpparams.runIdx    = 1;              
gpparams.datFormat = 'spikes';       % 'seq' if already binned and sqrt transf data
gpparams.method    = 'gpfa';         % 'gpfa' 'fa 'ppca' 'pca'
gpparams.xDim      = [10];             % Select number of latent dimensions
% JEB7 4-29 best xDim = 10
gpparams.binWidth  = 5;             % spike bin width in ms 
gpparams.numFolds  = 0;              % cross validation folds to find best xDim

% Extract neural trajectories
result = neuralTraj(gpparams.runIdx, dat,'datFormat', gpparams.datFormat, ...
                    'method', gpparams.method, 'xDims', gpparams.xDim, ...
                    'binWidth', gpparams.binWidth, 'numFolds', gpparams.numFolds);

% method = 'fa';
% result = neuralTraj(runIdx, dat,'datFormat', datFormat, ...
%                     'method', method, 'xDims', xDim, ...
%                     'binWidth', binWidth, 'numFolds', numFolds);

%%

% Orthonormalize neural trajectories
[estParams, seqTrain] = postprocess(result,'kernSD',80);
% NOTE: The importance of orthnormalization is described on 
%       pp.621-622 of Yu et al., J Neurophysiol, 2009.

gpfalatents = zeros(seqTrain(1).T,gpparams.xDim,numel(seqTrain));
ct = 1;
for i = numel(seqTrain):-1:1
    gpfalatents(:,:,ct) = seqTrain(1).xorth';
    ct = ct + 1;
end

%% save data

savepth = 'M:\Economo-Lab\data\gpfa';
savefn = [meta.anm '_' meta.date '_run' num2str(gpparams.runIdx) '.mat'];
save(fullfile(savepth,savefn),'gpfalatents','obj','params','meta','gpparams','-v7.3')

%% 

kernSD = 30; % select kernSD for two-stage methods
plotPredErrorVsDim(runIdx, kernSD);


%% plot projections

gpfatime = ((1:seqTrain(1).T).*(result.binWidth/1000)- mode(obj.bp.ev.(params.alignEvent)))./2;

f = figure;
for j = 1:gpparams.xDim
    ax = nexttile; hold on;
    for i = 1:numel(trials)

            temp = seqTrain(i);

            if ismember(trials(i), params.trialid{2})
                col = 'b';
            elseif ismember(trials(i), params.trialid{3})
                col = 'r';
            else
                continue
            end

            plot(gpfatime,mySmooth(temp.xorth(j,:),1),col)

    end
    hold off;
end


rhits = zeros(gpparams.xDim,seqTrain(1).T,numel(params.trialid{2}));
lhits = zeros(gpparams.xDim,seqTrain(1).T,numel(params.trialid{3}));
rix = 1;
lix = 1;
for i = 1:numel(trials)

    temp = seqTrain(i);

    if ismember(temp.trialId, params.trialid{2})
        rhits(:,:,rix) = temp.xorth;
        rix = rix + 1;
    elseif ismember(temp.trialId, params.trialid{3})
        lhits(:,:,rix) = temp.xorth;
        lix = lix + 1;
    else
        continue
    end
end

f = figure; hold on;
for j = 1:gpparams.xDim
    ax = nexttile; hold on;

    plot(gpfatime,mean(rhits(j,:,:),3),'b')
    plot(gpfatime,mean(lhits(j,:,:),3),'r')
    hold off;
end

 j%% Plot neural trajectories in 3D space

dims2plot = [1 5 6];

f = figure; hold on;
for i = 1:30%numel(trials)

    temp = seqTrain(i);

    if ismember(trials(i), params.trialid{2})
        col = 'b';
    elseif ismember(trials(i), params.trialid{3})
        col = 'r';
    else
        continue
    end

    plot3(temp.xorth(dims2plot(1),:),temp.xorth(dims2plot(2),:),temp.xorth(dims2plot(3),:),col)

end
hold off;


















