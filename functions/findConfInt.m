% Find confidence intervals for CDs Early, Late, and Go across sessions
% Find confidence intervals for ME across sessions
function [upperci,lowerci] = findConfInt(meta,taxis,conf,conditions)
nSessions = length(meta);
upperci.early = NaN(length(taxis)-1,2); lowerci.early = NaN(length(taxis)-1,2);
upperci.late = NaN(length(taxis)-1,2); lowerci.late = NaN(length(taxis)-1,2);
upperci.go = NaN(length(taxis)-1,2); lowerci.go = NaN(length(taxis)-1,2);
upperci.move = NaN(length(taxis)-1,2); lowerci.move = NaN(length(taxis)-1,2);
for c = 1:numel(conditions)
    upperci.early(:,c) = conf.cdearly.avg(2:end,c)+1.96*(conf.cdearly.std(2:end,c)/nSessions);  % Find the upper 95% confidence interval for each condition
    lowerci.early(:,c) = conf.cdearly.avg(2:end,c)-1.96*(conf.cdearly.std(2:end,c)/nSessions);  % Find lower 95% condifence interval for each condition

    upperci.late(:,c) = conf.cdlate.avg(2:end,c)+1.96*(conf.cdlate.std(2:end,c)/nSessions);  
    lowerci.late(:,c) = conf.cdlate.avg(2:end,c)-1.96*(conf.cdlate.std(2:end,c)/nSessions);  

    upperci.go(:,c) = conf.cdgo.avg(2:end,c)+1.96*(conf.cdgo.std(2:end,c)/nSessions);  
    lowerci.go(:,c) = conf.cdgo.avg(2:end,c)-1.96*(conf.cdgo.std(2:end,c)/nSessions);  

    upperci.move(:,c) = conf.movement.avg(2:end,c)+1.96*(conf.movement.std(2:end,c)/nSessions);  
    lowerci.move(:,c) = conf.movement.avg(2:end,c)-1.96*(conf.movement.std(2:end,c)/nSessions);  
end
end