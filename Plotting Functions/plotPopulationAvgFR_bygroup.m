% INPUTS:
% Ngroups = the number of groups that you want to split trials up into
% clr = colormap with Ngroups number of colors
% group = (1 x nTrials) array that assigns each trial a group
% ix = the sorted indices for how to group trials 
% groupsToPlot = which of the Ngroups you want to plot
% taxis = time axis
% groupslegend = the lengend 

% OUTPUTS:
% Will plot the population averaged firing rate for each group of trials
% (across cells and across all of the trials in the specified group)
function plotPopulationAvgFR_bygroup(obj,Ngroups,clr, group, ix,groupsToPlot,taxis,groupslegend)

for i = 1:Ngroups                 % For all groups of trials...
    c = clr(i, :);                                   % Find the color associated with the current group
    trix = ix(find(group==i));
    if ismember(i, groupsToPlot)                     % If current group is one that you want to plot...
        across_trials = mean(obj.trialpsth(:,:,trix),3); % Get mean across trials from this group
        ts = mean(medfilt1(across_trials, 25), 2,'omitnan');        % Find average FR from this group across cells
        plot(taxis,ts, 'Color', c, 'LineWidth', 3);
    end
end
xlim([-2.5 3])
legend(groupslegend,'Location','best')
xlabel('Time since go-cue (s)')
ylabel('Population avg FR (Hz)')
title('Population FR for each group')
hold off;

end % plotPopulationAvgFR_bygroup