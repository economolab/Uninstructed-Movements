function [upperci, lowerci] = getConfInt_Context(meta, avgCD, stdCD)
nSessions = length(meta);

upperci.AFC.true = avgCD.AFChit.true(2:end)+1.96*(stdCD.AFChit.true(2:end)/nSessions);  % Find the upper 95% confidence interval for each condition
lowerci.AFC.true = avgCD.AFChit.true(2:end)-1.96*(stdCD.AFChit.true(2:end)/nSessions);  % Find lower 95% condifence interval for each condition
upperci.FW.true = avgCD.FWhit.true(2:end)+1.96*(stdCD.FWhit.true(2:end)/nSessions);  
lowerci.FW.true = avgCD.FWhit.true(2:end)-1.96*(stdCD.FWhit.true(2:end)/nSessions);  

upperci.AFC.pred = avgCD.AFChit.pred(2:end)+1.96*(stdCD.AFChit.pred(2:end)/nSessions);  
lowerci.AFC.pred = avgCD.AFChit.pred(2:end)-1.96*(stdCD.AFChit.pred(2:end)/nSessions);  
upperci.FW.pred = avgCD.FWhit.pred(2:end)+1.96*(stdCD.FWhit.pred(2:end)/nSessions);  
lowerci.FW.pred = avgCD.FWhit.pred(2:end)-1.96*(stdCD.FWhit.pred(2:end)/nSessions);  
end