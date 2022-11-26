function [mode,dat,proj,kinfns] = calcKinModes(kin,obj,params,psthForProj,conditions,taxis)
kinfns = fieldnames(kin);
for i = 1:numel(kinfns)
    Y = kin.(kinfns{i});                                                    % feature data to use to calculate mode
    if ~strcmp(kinfns{i},'tongueVel') && ~strcmp(kinfns{i},'tongueAngle')
        e1 = find(taxis>-0.5,1,'first');
        e2 = find(taxis>-0.05,1,'first');
        params.tix = e1:e2;        % time points to use when finding mode (LATE DELAY)
    elseif strcmp(kinfns{i},'tongueVel') || strcmp(kinfns{i},'tongueAngle')
        e1 = find(taxis>0.05,1,'first');
        e2 = find(taxis>0.5,1,'first');
        params.tix = e1:e2;        % time points to use when finding mode (RESPONSE PERIOD)
    end

    if strcmp(kinfns{i},'tongueAngle') 
        conditions = 1:2;
    end
    [mode.(kinfns{i}), dat.(kinfns{i})] = findMode(obj, Y, params,conditions);         % Find the kinematic mode
    proj.(kinfns{i}) = getProjection(psthForProj, mode.(kinfns{i}));                   % Get the projection onto that mode
end
end  % calcKinModes