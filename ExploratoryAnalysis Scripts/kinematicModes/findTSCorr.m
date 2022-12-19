% TRIALS ARE CONCATENATED IN TIME
% kinmeas are observations (N x 1) -- kinematic measure over time (time x 1)
% PSTH are predictors (N x M)   -- activity of each cell over time (time x cells)
% r is loadings (M x 1)      -- corr coef of each cell with the kinematic feature

function tscorr = findTSCorr(PSTH, kinmeas,taxis,window,sortrange)
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
    FRepoch = permute(filtfr(e1:e2,:, :), [ 1 3 2]);                               % Re-order the dimensions of filtfr (switch the 2nd and 3rd dimensions)
    % And take only time-points that you want
    FRepoch = reshape(FRepoch, size(FRepoch, 1)*size(FRepoch, 2), size(FRepoch, 3));   % Concatenate trials in time (time*trials x cells)

    y = kinmeas(e1:e2, :);                                                  % Take time-points that you want of the kinematic feature
    y = reshape(y, [size(y, 1)*size(y, 2),1]);
    ix = find(isnan(y));
    y(ix) = 0;

    r = zeros(nCells, 1);       % (cells x 1)

    for i = 1:nCells            % For every cell...
        x = FRepoch(:, i);                       % Get the PSTH for all specified times and trials for that cell
        tmp = corrcoef(x, y);                    % Find the correlation coefficient between the cell's firing rate and the feature measure
        r(i) = tmp(1,2);                         % Store the corr coefficient (for this time range)
    end
    tscorr(:,count) = r;
end
delcorr = mean(tscorr(:,sortrange),2); [~,sortix] = sort(delcorr,'descend');
tscorr = tscorr(sortix,:);
end