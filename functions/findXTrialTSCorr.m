% TRIALS ARE CONCATENATED IN TIME
% kinmeas are observations (N x 1) -- kinematic measure over time (time x 1)
% PSTH are predictors (N x M)   -- activity of each cell over time (time x cells)
% r is loadings (M x 1)      -- corr coef of each cell with the kinematic feature

function tscorr = findXTrialTSCorr(PSTH, kinmeas,taxis,window,sortrange)
dt = abs(taxis(1)-taxis(2));    % Time jump between two time points
nSteps = round(window/dt);      % How many steps in the taxis you want to look over
plusmin = round(nSteps/2);      % How many steps before and after the current t you want to go
nCells = size(PSTH,2);
nWindows = length(((4+plusmin):(length(taxis) - plusmin)));
tscorr = NaN(size(PSTH,2),nWindows);            % (cells x number of time windows that you used)

fr = PSTH;                                      % Concatenate the single-trial PSTHs from the specified conditions (time x cells x trials)
fs = 1./mean(diff(taxis));                      % Find the sampling frequency for the time vector
fcut = 50;                                      % smoothing cutoff frequency

% Define and implement a Butterworth filter for the single-trial PSTHs
[b,a] = butter(2,fcut./fs./2);                  % Butterworth filter: 2nd order, cutoff frequency determined using fcut and the sampling rate
% b = filter coefficient vector of numerator
% a = filter coefficient vector of denominators

filtfr = filtfilt(b, a, fr);                    % Digital filtering of the single-trial PSTHs with the filter described by vectors A and B

count = 0;
for t = 4+plusmin:(length(taxis) - plusmin)                  % Start at the 5th index in the time axis because first few are usually NaNs
    count = count+1;
    e1 = t-plusmin;
    e2 = t+plusmin;
    FRepoch = filtfr(e1:e2,:, :);                                  % Get the firing rates during the epoch that you want
    meanFR = squeeze(mean(FRepoch,1,'omitnan'));                   % Get the avg FR for all neurons during this epoch (on each trials)
    
    y = kinmeas(e1:e2, :);                                         % Take time-points that you want of the kinematic feature
    ix = find(isnan(y));
    y(ix) = 0;
    meanKin = mean(y,1,'omitnan');                                 % Get the avg kin measure during this epoch on each trial 
    
    r = NaN(nCells,1);
    for i = 1:nCells            % For every cell...
        currCellFR = meanFR(i,:);
        tmp = corrcoef(currCellFR, meanKin);                       % Find the correlation coefficient between the cell's avg firing rate
                                                                   % during this epoch and the avg feature measure 
        r(i) = tmp(1,2);
    end
    tscorr(:,count) = r;
end
delcorr = mean(tscorr(:,sortrange),2); [~,sortix] = sort(delcorr,'descend');
tscorr = tscorr(sortix,:);
end