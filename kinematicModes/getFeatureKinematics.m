function [featPos,featSpeed] = getFeatureKinematics(taxis,obj,conditions,met,view,feat,params)

[xpos, ypos] = findPosition(taxis, obj, conditions, met, view, feat,params);
[xvel, yvel] = findVelocity(taxis, obj, conditions, met, view, feat,params);
featPos = [];
featSpeed = [];
for i = 1:numel(conditions)
    featPos = [featPos sqrt(xpos{i}.^2+ypos{i}.^2)];            % Find displacement 
    featSpeed = [featSpeed sqrt(xvel{i}.^2+yvel{i}.^2)];        % Find speed in x and y direction
end