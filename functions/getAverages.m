% Function for finding the average value of a measurement within a specific
% time interval

% INPUTS: timeInt = the time interval (in indices, within the trial) that you want to average
% over
% allTrials = [time x nTrials] array where you have values for all points
% in time across a trial, for all trials 

% OUTPUT: [1 x nTrials] array where you have a single averaged value of
% your measurement for each trial
function avgValue_acrossTrials = getAverages(timeInt,allTrials)

% Find average jaw velocity during given time interval for each trial 
avgValue_acrossTrials = mean(allTrials(timeInt, :), 1); 

end %getAverages