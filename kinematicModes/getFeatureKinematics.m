function [featPos,featSpeed] = getFeatureKinematics(taxis,obj,conditions,met,view,feat)

[xpos, ypos] = findPosition(taxis, obj, conditions, met, view, feat);
[xvel, yvel] = findVelocity(taxis, obj, conditions, met, view, feat);
featPos = [sqrt(xpos{1}.^2+ypos{1}.^2) sqrt(xpos{2}.^2+ypos{2}.^2)];
featSpeed = [sqrt(xpos{1}.^2+ypos{1}.^2) sqrt(xpos{2}.^2+ypos{2}.^2)];

end