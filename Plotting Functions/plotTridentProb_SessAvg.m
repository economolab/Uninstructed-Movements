% Plots the probability of jaw movement during a session, separated by
% trial type

% INPUTS: obj, met for the current session
% conditions: cell array (1 x num conditions) the behavioral conditions that you want to look at jaw
% movement for 
% colors: cell array containing the colors that you want to be used in the
% plot 
function plotTridentProb_SessAvg(obj,met,conditions,colors)
edges = met.tmin:met.dt:met.tmax;

tridentprob = tridentProbSessionAvg(obj,met,conditions);

for i = 1:numel(tridentprob)
    plot(edges, tridentprob{i},'color', colors{i}, 'LineWidth', 2);
    hold on;
end

xlabel('Time since go-cue (s)','FontSize',13)
ylabel('Prob of jaw movement','FontSize',13)
xlim([-2.3 2.5])

end  % end 'plotJawProb_SessAvg'