% Function for finding the average probability of movement at all time
% points in a trial

% INPUTS: 
% conditions = which trial conditions you want to look at
% mov = whether animal is moving at any points during trial
% me = raw motion energy values for whole session

% OUTPUT: moveprob: cell array (1 x num conditions)
% Each cell will contain the probability of movement for one specified
% condition

function moveprob = moveProbSessionAvg(met,conditions,mov,me,obj,params,taxis)

moveprob = cell(1,numel(conditions));
for i = 1:numel(conditions)
    cond = conditions{i};
    trialsToUse = met.trialid{cond};

    moveThresh = me.moveThresh;  % If the motion energy crosses this threshold, the animal is considered to be moving
    edges = taxis;
    Ntrials = length(trialsToUse);

    move = NaN(numel(edges), length(trialsToUse));

    for ii = 1:length(trialsToUse)
        q = trialsToUse(ii);
        trialtime = [mov.moveTime{q};mov.stationaryTime{q}];        % Get all the times in the current trial
        trialtime = sort(trialtime,'ascend'); trialtime = trialtime - obj.bp.ev.(params.alignEvent)(q);
        
        ts = me.data{q};                                            % Motion energy values for this trial
        tsinterp = interp1(trialtime,ts, edges);                    % Linear interpolation of ME to keep number of time points consistent across trials

        %Find when the difference between the jaw velocity and the
        %baseline jaw
        %velocity is above a given threshold (when is jaw moving?)
        move(:, ii) = tsinterp;%>moveThresh;                      % Store binary array of whether or not animal is moving
    end
    
    
    movedat = mean(move,2,'omitnan');         % Take the average probability of jaw movement across trials
    movedat = mySmooth(movedat,21);           % Smooth the avg probability of jaw movement 
    moveprob{i} = movedat;                    % Store the avg prob of jaw movement for the current condition

end
%end % moveProbSessionAvg