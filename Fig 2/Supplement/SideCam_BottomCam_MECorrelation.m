% Compare motion energy found via the side cam vs bottom cam
clear,clc,close all

whichcomp = 'LabPC';                                                % LabPC or Laptop

% Base path for code depending on laptop or lab PC
if strcmp(whichcomp,'LabPC')
    basepth = 'C:\Users\Jackie Birnbaum\Documents\Code';
elseif strcmp(whichcomp,'Laptop')
    basepth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab\Code';
end

% add paths for data loading scripts, all fig funcs, and utils
utilspth = [basepth '\Munib Uninstruct Move\uninstructedMovements_v2'];
addpath(genpath(fullfile(utilspth,'DataLoadingScripts')));
addpath(genpath(fullfile(utilspth,'funcs')));
addpath(genpath(fullfile(utilspth,'utils')));
addpath(genpath(fullfile(utilspth,'fig1')));
figpth = [basepth '\Uninstructed-Movements\Fig 2'];
addpath(genpath(fullfile(figpth,'funcs')));
addpath(genpath(fullfile(figpth,'Hazarded Delay')));
%% LOAD Hazarded Delay data
% Everything aligned to the delay period; tmin = -2.5, tmax = 3; dt = 1/200
% HAVE TO ACCOUNT FOR 0.5 SECOND DELAY WITHOUT SPIKEGLX
load([basepth '\Uninstructed-Movements\Fig 2\Hazarded Delay\BehavData\JEB11_JEB12_objs.mat']);
load([basepth '\Uninstructed-Movements\Fig 2\Hazarded Delay\BehavData\JEB11_JEB12_meta.mat']);
load([basepth '\Uninstructed-Movements\Fig 2\Hazarded Delay\BehavData\JEB11_JEB12_params.mat']);

% Set params for hazarded delay data
params.tmin = meta(1).tmin; params.tmax = meta(1).tmax;
params.traj_features = {{'tongue','left_tongue','right_tongue','jaw','trident','nose'},...
    {'top_tongue','topleft_tongue','bottom_tongue','bottomleft_tongue','jaw'}};
params.smooth = 15;
params.advance_movement = 0;
params.feat_varToExplain = 80; % num factors for dim reduction of video features should explain this mu

% Re-structure Hazarded Delay data to be the same format
[params,obj] = cleanUpData(params, meta,objs);

%% Load motion energy for hazarded delay sessions
if strcmp(whichcomp,'LabPC')
    datapth = 'C:\Users\Jackie Birnbaum\Documents\Data';
elseif strcmp(whichcomp,'Laptop')
    datapth = 'C:\Users\Jackie\Documents\Grad School\Economo Lab';
end

% Side cam
camview = 'sidecam';
for sessix = 1:numel(meta)
    disp(['Loading side cam ME for session ' num2str(sessix)])
    if obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{1},2) || obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{2},2)
        disp('Number of Bpod trials does not match number of video files for this session')
    else
        sideme(sessix) = loadMotionEnergy_Behav(obj(sessix), meta(sessix), params(sessix), datapth,camview);
    end
end

% Bottom cam
camview = 'bottomcam';
for sessix = 1:numel(meta)
    disp(['Loading bottom cam ME for session ' num2str(sessix)])
    if obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{1},2) || obj(sessix).bp.Ntrials ~= size(obj(sessix).traj{2},2)
        disp('Number of Bpod trials does not match number of video files for this session')
    else
        bottomme(sessix) = loadMotionEnergy_Behav(obj(sessix), meta(sessix), params(sessix), datapth,camview);
    end
end
%% Method 1: On each trial, find the R^2 between the side cam ME and the bottom cam ME
% Then take the average across all trials to find one R^2 value for each session
sm = 51;
M1_R2 = NaN(1,length(meta));
for sessix = 1:length(meta)
    % Check to make sure side cam and bottom cam have same number of trials
    if size(sideme(sessix).data,2)==size(bottomme(sessix).data,2)
        nTrials = size(sideme(sessix).data,2);
        tempR2 = NaN(1,nTrials);
        for trix = 1:nTrials
            side = mySmooth(sideme(sessix).data(:,trix),sm);
            bottom = mySmooth(bottomme(sessix).data(:,trix),sm);
            R2 = corrcoef(side,bottom);
            tempR2(trix) = R2(2);
%             plot(side); hold on; plot(bottom); hold off;
%             title(num2str(tempR2(trix)))
%             pause
        end
        M1_R2(sessix) = mean(tempR2,'omitnan');
    end
end
%% Method 2: Concatenate the side cam ME from all trials together into one long time series
% Do the same for bottom cam ME.  Find the correlation between those two long time series 
M2_R2 = NaN(1,length(meta));
for sessix = 1:length(meta)
    % Check to make sure side cam and bottom cam have same number of trials
    if size(sideme(sessix).data,2)==size(bottomme(sessix).data,2)
        nTrials = size(sideme(sessix).data,2);
        side_alltrix = [];
        bottom_alltrix = [];
        for trix = 1:nTrials
            side = mySmooth(sideme(sessix).data(:,trix),sm);
            side_alltrix = [side; side_alltrix];
            bottom = mySmooth(bottomme(sessix).data(:,trix),sm);
            bottom_alltrix = [bottom; bottom_alltrix];  
        end
        R2 = corrcoef(side_alltrix,bottom_alltrix);
%         plot(side_alltrix); hold on; plot(bottom_alltrix); hold off;
%         pause
        M2_R2(sessix) = R2(2);
    end
end
%% Plot histograms of R2 values found using either method
nBins = 5;
subplot(2,1,1)
histogram(M1_R2,nBins)
title('Method 1: R2 on each trial, avg all')
xlabel('R2 val')
ylabel('Num sessions')

subplot(2,1,2)
histogram(M2_R2,nBins)
title('Method 2: Concatenate all trials, R2 once')
xlabel('R2 val')
ylabel('Num sessions')
%% Compare percentage of time spent moving found via each camera
sess = [1,3:6,8:14];
x = linspace(0, 1, 100);
y = linspace(0, 1, 100);
for sessix = sess
    % Check to make sure side cam and bottom cam have same number of trials
    if size(sideme(sessix).data,2)==size(bottomme(sessix).data,2)
        nTrials = size(sideme(sessix).data,2);
        nTimepts = size(sideme(sessix).data,1);
        sidepct = NaN(nTrials, 1);
        bottompct = NaN(nTrials, 1);
        for trix = 1:nTrials
            side = sum(sideme(sessix).data(:,trix)>sideme(sessix).moveThresh);       % Get number of time points where sidecam ME is above movethreshold
            sidepct(trix) = side/nTimepts;                                           % Divide by total number of time points to get percentage
            bottom = sum(bottomme(sessix).data(:,trix)>bottomme(sessix).moveThresh); % Get number of time points where bottomcam ME is above movethreshold
            bottompct(trix) = bottom/nTimepts;
        end
        nexttile
        scatter(sidepct,bottompct,10,'filled'); hold on;
        plot(x,y,'k--')
        xlim([0 1])
        xlabel('side cam')
        ylim([0 1])
        ylabel('bottom cam')
        title(['Pct of time moving; Session' num2str(sessix)])
    end
end

%% FUNCTIONS
function [params,obj] = cleanUpData(params, meta,objs)
% Reorganize params into a struct to match structure of updated data
temp = params;
clear params
for sessix = 1:numel(meta)
    temp2 = temp;
    temp2.trialid = meta(sessix).trialid;
    params(sessix) = temp2;
end
clear temp temp2

% Reorganize objs into a struct to match structure of updated data
temp = objs;
clear objs
for sessix = 1:numel(meta)
    temp2 = temp{sessix};
    temp2.time = params(sessix).taxis; %+0.5
    obj(sessix) = temp2;
end
end