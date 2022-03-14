clear,clc,close all

% add folders in current directory to the path
% this way you have access to functions in funcs and utils directories
addpath(genpath(pwd))


%% SET PARAMS
params.alignEvent          = 'jawOnset'; % 'goCue'  'jawOnset'  'firstLick'  'lastLick'

% time warping only operates on neural data for now.
% TODO: time warp for video and bpod data
params.timeWarp            = 0;  % piecewise linear time warping - each lick duration on each trial gets warped to median lick duration for that lick across trials
params.nLicks              = 20; % number of post go cue licks to calculate median lick duration for and warp individual trials to 

params.lowFR               = 1; % remove clusters with firing rates across all trials less than this val

% set conditions to calculate PSTHs for
params.condition(1)     = {'R&hit&~stim.enable&~autowater&~early'};         % right hits, no stim, aw off
params.condition(end+1) = {'L&hit&~stim.enable&~autowater&~early'};         % left hits, no stim, aw off
params.condition(end+1) = {'R&hit&~stim.enable&autowater&~early'};          % right hits, no stim, aw on
params.condition(end+1) = {'L&hit&~stim.enable&autowater&~early'};          % left hits, no stim, aw on
% params.condition(end+1) = {'R&miss&~stim.enable&~autowater&~early'};        % error right, no stim, aw off
% params.condition(end+1) = {'L&miss&~stim.enable&~autowater&~early'};        % error left, no stim, aw off
% params.condition(end+1) = {'~hit&~miss&~stim.enable&~autowater&~early'};    % ignore
% params.condition(end+1) = {'hit&~stim.enable&~autowater&~early'};           % hit 2afc
% params.condition(end+1) = {'hit&~stim.enable&autowater&~early'};            % hit aw

% specify probe number and areas to load and process data
params.probe(1) = 1;
params.probeArea{1} = 'ALM';
% 
% params.probe(end+1) = 2;
% params.probeArea{end+1} = 'ALM';

params.tmin = -2.5;
params.tmax = 2.5;
params.dt = 1/200;

% smooth with causal gaussian kernel
params.smooth = 15;

% cluster qualities to use
params.quality = {'all'}; % accepts any cell array of strings - special character 'all' returns clusters of any quality

% video Hz
params.vid_dt = 1/400;

%% SET METADATA
% experiment meta data
meta.datapth = '/Users/Munib/Documents/Economo-Lab/code/data';
% meta.datapth = '/Volumes/MUNIB_SSD/Experiments';
meta.anm = 'JGR2';
meta.date = '2021-11-16';
meta.datafn = findDataFn(meta);

%% LOAD AND PROCESS DATA
dat = load(fullfile(meta.datapth,meta.datafn));
obj = dat.obj;

for prbnum = 1:numel(params.probe)
    disp('______________________________________________________')
    disp(['Processing data for probe ' num2str(prbnum)])
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

%% find lick times 
if ~isfield(obj,'lick') % this is true when params.timeWarp = 0
    % get first params.Nlicks lick times, for each trial
    [obj.lick.lickStart,obj.lick.lickEnd,obj.lick.lickDur] = findLickTimes(obj,params.nLicks);
    %find median lick times for each lick across trials
    % only using trials where a lick was present in median calculation
    obj.lick.med = findMedianLickTimes(obj.lick.lickStart,obj.lick.lickEnd,obj.lick.lickDur, params.nLicks);
end

%% find time to nth post go cue lick
% 
% plt.legend = {'R 2AFC','L 2AFC','R AW', 'L AW'};
% plt.colors = {[0 0 1],[1 0 0], ...
%                  [190, 3, 252]./255,[252, 190, 3]./255};
% plt.save = 1;
%              
% conds = [1 2 3 4];
% lickNums = [1:10];
% for i = 1:numel(lickNums)
%     lickNum = lickNums(i);
%     lt = findLickTimesByCond(conds,params,obj.lick.lickStart,lickNum, plt);
% end


%% analyze jaw traces for each trial type from bottom cam
% analyzeJaw(obj);

%% calculate eigenpose
% calcEigenpose(obj);

%% get behavioral trajectories

featNames = getFeatNames(obj);

feat2use = {[4 5 6],
             []}; % which features to use for each view. If empty array, won't use that view

% feature trajectories
traj = getFeatTraj(obj,featNames,feat2use);

% align feature trajectories to move onset (estimated with jaw position)
% trials that don't have t_preEv or t_postEv time points to use are
% discarded. See trialid field for the actual trial number. 
t_preEv = 2;
t_postEv = 2;
traj = alignFeatTraj(traj,obj.bp.ev.jawOnset,t_preEv,t_postEv,params.vid_dt);

% feature velocities
sm = 25; % smoothing window for sgolay filter
vel = getFeatVelocity(traj,params.vid_dt,sm);

% % plot velocity for single trial
% for trix = 1:numel(traj)
%     plotFeatVelocity(traj,vel,trix);
%     pause
% end

%% analysis

% use a single trial for now
trix = 99;

% get data for trial into a (time,features) matrix
time = traj(trix).time;

[datmat,labels] = featStruct2Mat(vel(trix));

%% k-means clustering
%  https://sleap.ai/notebooks/Analysis_examples.html

nstates = 3;

idx = kmeans(datmat,nstates);

plotKMeansEthogram(trix,time,datmat,labels,idx);


%% more sophisiticated clustering

% procedure:
% 1. wavelet transform w/ morlet wavelet (time-frequency decomp). 
% 2. tsne 2d embedding
% 3. watershed segmentation (clustering)

% refs:
% 1. https://royalsocietypublishing.org/doi/10.1098/rsif.2014.0672
% 2. https://www.nature.com/articles/s41592-018-0234-5#Sec8
% 3. https://www.sciencedirect.com/science/article/pii/S0896627320308941#sec5


%% morlet wavelet transform

%%%%%%%% Wavelet Parameters %%%%%%%%
%number of wavelet frequencies to use
numPeriods = 25;
%dimensionless Morlet wavelet parameter
omega0 = 5;
%minimum frequency for wavelet transform (Hz)
minF = 0.5;
%maximum frequency for wavelet transform (Hz)
maxF = 60;

minT = 1 ./ maxF;
maxT = 1 ./ minF;
Ts = minT.*2.^((0:numPeriods-1).*log(maxT/minT)/(log(2)*(numPeriods-1)));
f = fliplr(1./Ts);


for featix = 1:size(datmat,2)
    amplitudes{featix} = fastWavelet_morlet_convolution_parallel(datmat(:,featix),f,omega0,params.vid_dt);
end

plotCWTAmplitude(time,amplitudes,f,labels,trix);

%% tsne 2d embedding

% initialize parameters
parameters = set_tSneParameters();

traindat = cell2mat(amplitudes'); % (nFeat*nFreq,time)
% imagesc(traindat)

% functionality
[trainingEmbedding,betas,P,errors] = run_tSne(traindat',parameters);


G = zeros(size(time));
G(ceil(numel(time)/2):end) = 1;
figure;
subplot(2,1,1)
scatter(Y(:,1),Y(:,2),20,time,'filled'); colorbar;
subplot(2,1,2)
scatter(Y(:,1),Y(:,2),20,G,'filled');


% kmeans with traindat
nstates = 4;

idx = kmeans(traindat,nstates);

plotKMeansEthogram(trix,time,datmat,labels,idx);






