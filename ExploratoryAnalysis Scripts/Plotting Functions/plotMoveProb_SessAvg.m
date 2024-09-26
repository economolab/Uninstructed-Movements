% Plots the probability of movement during a session, separated by
% trial type

% INPUTS: obj, met for the current session
% conditions: cell array (1 x num conditions) the behavioral conditions that you want to look at movement for 
% colors: cell array containing the colors that you want to be used in the
% plot 
function plotMoveProb_SessAvg(met,conditions,colors,mov,me,obj,params,taxis)
edges = taxis;

moveprob = moveProbSessionAvg(met,conditions,mov,me,obj,params,taxis);

for i = 1:numel(moveprob)
    plot(edges, moveprob{i},'color', colors{i}, 'LineWidth', 2);
    hold on;
end

end  % 'plotMoveProb_SessAvg'