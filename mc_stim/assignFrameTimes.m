% DON'T USE: just needed this temporarily
for i = 1:numel(obj)
    if isfield(obj(i).traj{1}(1),'frameTimes')
        if ~isfield(obj(i).traj{1}(1),'NdroppedFrames')
            for j = 1:obj(i).bp.Ntrials
                obj(i).traj{1}(j).NdroppedFrames = 0;
                obj(i).traj{2}(j).NdroppedFrames = 0;
            end
        end
        continue
    end
    for j = 1:obj(i).bp.Ntrials
        temp = obj(i).traj{1}(j);
        frameTimes = (0:(size(temp.ts,1)-1))./400;
        obj(i).traj{1}(j).frameTimes = frameTimes;
        obj(i).traj{2}(j).frameTimes = frameTimes;
        obj(i).traj{1}(j).NdroppedFrames = 0;
        obj(i).traj{2}(j).NdroppedFrames = 0;
    end
end