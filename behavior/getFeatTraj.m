function traj = getFeatTraj(obj,featNames,feat2use)

traj = struct();

% get trajectories
ct = 1;
for vix = 1:numel(feat2use) % for each view
    
    
    nTrials = numel(obj.traj{vix});
    for trix = 1:nTrials % for each trial in view
        traj(trix).trialid = trix;
        traj(trix).time = obj.traj{vix}(trix).frameTimes-0.5;
        for fix = 1:numel(feat2use{vix}) % for each feature to use for current view
            fname = featNames{vix}{feat2use{vix}(fix)};
            traj(trix).([fname '_cam' num2str(vix-1)]) = obj.traj{vix}(trix).ts(:,:,feat2use{vix}(fix));
        end
    end
    
    
end



end