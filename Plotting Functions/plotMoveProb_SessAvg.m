% Plots the probability of movement during a session, separated by
% trial type

% INPUTS: obj, met for the current session
% conditions: cell array (1 x num conditions) the behavioral conditions that you want to look at movement for 
% colors: cell array containing the colors that you want to be used in the
% plot 
function plotMoveProb_SessAvg(met,conditions,colors,mov,me)
edges = met.tmin:met.dt:met.tmax;

moveprob = moveProbSessionAvg(met,conditions,mov,me);

for i = 1:numel(moveprob)
    plot(edges, moveprob{i},'color', colors{i}, 'LineWidth', 2);
    hold on;
end

xlabel('Time since go-cue (s)','FontSize',13)
ylabel('Prob of movement','FontSize',13)
xlim([-2.3 2.5])

end  % 'plotMoveProb_SessAvg'