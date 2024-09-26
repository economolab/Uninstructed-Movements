function [mode,dat,proj,kinfns] = calcKinModes(kin,obj,params,psthForProj,conditions)

for feat = 1:numel(kinfns)
    Y = kin.(kinfns{feat});                                                    % feature data to use to calculate mode
    if ~strcmp(kinfns{feat},'tongueVel') && ~strcmp(kinfns{feat},'tongueAngle')
        e1 = find(obj.time>-0.5,1,'first');
        e2 = find(obj.time>-0.05,1,'first');
        params.tix = e1:e2;        % time points to use when finding mode (LATE DELAY)
    elseif strcmp(kinfns{feat},'tongueVel') || strcmp(kinfns{feat},'tongueAngle')
        e1 = find(obj.time>0.05,1,'first');
        e2 = find(obj.time>0.5,1,'first');
        params.tix = e1:e2;        % time points to use when finding mode (RESPONSE PERIOD)
    end

    if strcmp(kinfns{feat},'tongueAngle') 
        conditions = 1:2;
    end
    [mode.(kinfns{feat}), dat.(kinfns{feat})] = findMode(obj, Y, params,conditions);         % Find the kinematic mode
    proj.(kinfns{feat}) = getProjection(psthForProj, mode.(kinfns{feat}));                   % Get the projection onto that mode
end
end  % calcKinModes