function [featPos,featSpeed] = getTongueKinematics(taxis,obj,conditions,met,view,feat)

[xpos, ypos] = findTonguePosition(taxis, obj, conditions, met, view, feat);
[xvel, yvel] = findTongueVelocity(taxis, obj, conditions, met, view, feat);
featPos = [];
featSpeed = [];
for i = 1:numel(conditions)
    featPos = [featPos sqrt(xpos{i}.^2+ypos{i}.^2)];            % Find displacement 
    featSpeed = [featSpeed sqrt(xvel{i}.^2+yvel{i}.^2)];        % Find speed in x and y direction
end