% Plots the probability of jaw movement during a session, separated by
% trial type

% INPUTS: obj, met for the current session
% conditions: cell array (1 x num conditions) the behavioral conditions that you want to look at jaw
% movement for 
% colors: cell array containing the colors that you want to be used in the
% plot 
function plotJawVel_SessAvg(obj,met,conditions,colors,taxis,conf)
edges = met.tmin:met.dt:met.tmax;

[jawprob,jawstd] = jawVelSessionAvg(obj,met,conditions);

nCond = numel(jawprob);
upperci = NaN(length(taxis),nCond);
lowerci = NaN(length(taxis),nCond);
for c = 1:numel(jawprob)
    temp = jawprob{c};
    tempstd = jawstd{c};
    plot(edges, temp,'color', colors{c}, 'LineWidth', 3);
    hold on;
    nTrials = size(temp,2);
    upperci(:,c) = temp(2:end)+1.96*(tempstd(2:end)/nTrials);  % Find the upper 95% confidence interval for each condition      
    lowerci(:,c) = temp(2:end)-1.96*(tempstd(2:end)/nTrials);  % Find lower 95% condifence interval for each condition
    if strcmp(conf,'yes')
        patch([taxis(9:end) fliplr(taxis(9:end))],[lowerci(9:end,c)' fliplr(upperci(9:end,c)')],colors{c},'FaceAlpha',0.2,'EdgeColor','none')
    end
end

xlabel('Time since go-cue (s)','FontSize',13)
ylabel('Jaw velocity','FontSize',13)
xlim([-2.3 2.5])

end  % end 'plotJawProb_SessAvg'