% Function for identifying the indices of trials where the animal is moving
% its jaw 15% of the delay period (i.e. Early Jaw trial)

% Jaw is moving if the variance of jaw velocity crosses a given threshold

function obj = findEarlyJaw(obj)   
currtraj = obj.traj{1};         % Get side-view video trajectories from current data object
earlyJawix = [];
thresh = 0.12;                  % Threshold for variance of jaw velocities
for j = 1:size(currtraj,2)           % Go through all trials...
    Jawx = mySmooth(currtraj(j).ts(:,1,2),9);    % Get x-values for jaw trajectories for current trial
    ts = diff(mySmooth(Jawx, 21));               % Get jaw velocity and smooth it
    sig = movvar(ts, 100);                       % Find the variance of the jaw velocity at all time points (using a sliding window)
    sig = [0; sig];
    goCue = obj.bp.ev.goCue(j);              % Find time of goCue for current trial
    sample = obj.bp.ev.sample(j)+0.05;       % Find time of sample for current trial
    t = currtraj(j).frameTimes - 0.5;            % Get trial time for all frames
    %figure(i+numel(objs));
    %         plot(t(1:end-1),3*j+ts); hold on; plot(t,3*j+sig);
    
    ix = find(sig>thresh & t<goCue & t>sample);  % Find indices during delay where the variance is greater than the threshold
    
    delayix = find(t<(goCue-0.05) & t>sample);   % Find all indices during delay period 
    moveThresh = 0.15*(length(delayix));         % Find number of frames that = 15% of delay epoch frames

    if length(ix) > moveThresh                   % If the jaw is moving in more than 15% of delay epoch frames...
        earlyJawix = [earlyJawix,j];             % Store the trial number as an early move trial
    end
end
obj.earlyJawix = earlyJawix';             % Store the trial numbers for current object that have early jaw movements
end  %findEarlyJaw