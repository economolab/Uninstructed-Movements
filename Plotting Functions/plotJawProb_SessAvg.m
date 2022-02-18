% Plots the probability of jaw movement during a session, separated by
% trial type

% INPUTS: dirs = Cell array containing the trial type you want to plot (R or
% L)
% contexts = Cell array containing the context that you want to plot for
% ('afc' or 'auw' for 2AFC or Autowater)

function plotJawProb_SessAvg (dirs,contexts,anm,date,colors)
jawdat = cell(1,numel(dirs));
for ii = 1:numel(dirs)
    jawdat{ii} = jawkinSessionAvg(dirs{ii},contexts{ii},anm,date);
    y = mean(jawdat{ii},2,'omitnan');
    plot(edges, y(2:end),'color', colors{ii}, 'LineWidth', 2);
    hold on;
end
xlabel('Time since go-cue (s)','FontSize',13)
ylabel('Prob of jaw movement','FontSize',13)
xlim([-2.3 2.5])