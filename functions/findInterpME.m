% Find motion energy at all points during the trial

function [MEinterp,MEraw] = findInterpME(edges,conditions, met,mov,me)
MEinterp = cell(1,numel(conditions));
MEraw = cell(1,numel(conditions));

for c = 1:numel(conditions)
    cond = conditions{c};
    nTrials = numel(met.trialid{cond});
    tempmove = nan(numel(edges)-1, nTrials);    % (time x num trials in curr condition)
    tempme = nan(numel(edges)-1, nTrials);    % (time x num trials in curr condition)

    for i = 1:nTrials                        % For every trial in the condition
        q = met.trialid{cond}(i);
        
        trialtime = [mov.moveTime{q};mov.stationaryTime{q}];        % Get all the times in the current trial
        trialtime = sort(trialtime,'ascend');
        
        ts = me.data{q};                                            % Motion energy values for this trial
        tsinterp = interp1(trialtime,ts, edges);                    % Linear interpolation of ME to keep number of time points consistent across trials

        tempmove(:,i) = tsinterp(2:end);                                   % Store interpolated ME for the trial
        tempme(:,i) = ts(2:length(edges));                          % Store raw ME for the trial
    end

    MEinterp{c} = tempmove;
    MEraw{c} = tempme;
end

end  % findInterpME