function traj_aligned = alignFeatTraj(traj,jawOnset,tpre,tpost,dt)

traj_aligned = traj;
discard = zeros(size(traj));
for trix = 1:numel(traj) % for each trial
    
    fnames = fieldnames(traj(trix));
    
    % for current trial, find ix to keep
    onset = jawOnset(trix);
    [~,onsetix] = min(abs(traj(trix).time - onset));
    ix1 = onsetix - ceil(tpre/dt);
    ix2 = onsetix + ceil(tpost/dt);
    if ix1 < 1 || ix2 > numel(traj(trix).time)
        discard(trix) = 1;
        continue
    end
    
    for fix = 1:numel(fnames) % for each field in traj
        % trim each time series in the struct
        if strcmpi(fnames{fix},'trialid')
            continue
        elseif strcmpi(fnames{fix},'time')
            traj_aligned(trix).(fnames{fix}) = traj_aligned(trix).(fnames{fix})(ix1:ix2,:) - onset;
            continue
        end
        traj_aligned(trix).(fnames{fix}) = traj(trix).(fnames{fix})(ix1:ix2,:);
        
    end
    
end

traj_aligned = traj_aligned(~discard);


end