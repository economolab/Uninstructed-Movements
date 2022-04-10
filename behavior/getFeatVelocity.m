function vel = getFeatVelocity(traj,dt,sm)

vel = struct();

for trix = 1:numel(traj) % for each trial
    
    fnames = fieldnames(traj(trix));
    for fix = 3:numel(fnames) % for each field in traj that's not trialid or time
        
        % set vel to position
        vel(trix).(fnames{fix}) = traj(trix).(fnames{fix})(:,[1 2]);
        
        for cix = 1:2 % for each coordinate (x/y)
            % calculate velocity for feature and coordinate
            try
            vel(trix).(fnames{fix})(:,cix) = myDiff(vel(trix).(fnames{fix})(:,cix),dt);
            catch
                'a'
            end
            % smooth velocity sgolayfilt(x,order,framelen)
            vel(trix).(fnames{fix})(:,cix) = sgolayfilt(vel(trix).(fnames{fix})(:,cix),3,sm);
%             % get norm of velocity (this is what SLEAP does)
%             vel(trix).(fnames{fix})(:,cix) = norm(vel(trix).(fnames{fix})(:,cix));
        end
    end
    
end



end