% Plots the probability of jaw movement for right and left trials,
% separated by the length of the delay period

% INPUTS:
% taxis = (time x 1) vector where 0 is go-cue
% jawvel = struct array with two fields (left and right).  Each field is a
% (1 x numDelays) cell array.  Each cell contains (time x 1) vector of jaw
% vel
% params = params variable 
function plotJawProb_HazardDel(taxis, jawvel,params)
for g = 1:length(params.delay)                  % For the first four delay lengths...
    subplot(2,2,g)
    plot(taxis,jawvel.right{g},'LineWidth',1.5); hold on;   % Plot the avg jaw prob on the right for that delay length
    plot(taxis,jawvel.left{g},'LineWidth',1.5); hold on;    % Same for left
    ylim([0 1])
    xlim([-2.5 2])
    xline(0-params.delay(g)-1.3,'LineStyle','-.','LineWidth',1) % Draw vert line at sample start
    xline(0-params.delay(g),'LineStyle','-.','LineWidth',1)     % Draw vert line at sample end
    xline(0,'LineStyle','--','LineWidth',1)                     % Draw vert line at go-cue
    legend('Right','Left')
    title('Delay =',params.delay(g),'FontSize',13)
    xlabel('Time since go-cue (s)')
    ylabel('Probability of jaw movement (%)')
end
end   % end plotJawProb_HazardDel