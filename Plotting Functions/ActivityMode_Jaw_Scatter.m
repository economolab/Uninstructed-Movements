% Plots a scatter plot of jaw velocity vs. choice mode on individual trials
% Each data point = an average value of jaw vel or activity mode during a
% specific time epoch

% INPUTS: jawVel = [1 x nTrials] array of average jaw velocities
% avgChoice = [1 x nTrials] array of average choice mode values
% conditions = range of numbers to indicate which conditions you want to
% make this plot for (i.e. 1:2 will plot R and L 2AFC hits)
% clrs = cell array of colors to use

function ActivityMode_Jaw_Scatter(jawVel,avgChoice,conditions,met,clrs,obj,params)
cnt = 0;
for j = conditions                              % For all conditions                 
    trialid = met.trialid{j};                   % Get the trials that correspond to those conditions
    nTrials = length(trialid);
    if cnt == 0
        TrialsToUse = cnt+1:nTrials;
    else
        TrialsToUse = cnt:cnt+nTrials-1; 
    end
    cnt = cnt+1+nTrials;

    color = clrs{j};                     
    scatter(jawVel(TrialsToUse),avgChoice(TrialsToUse),'MarkerFaceColor',color)     % Make a scatter plot comparing these values                                                                       % for the desired trials, colored according to condition
    hold on;
end
ax = gca;
ax.FontSize = 12;
end  %ActivityMode_Jaw_Scatter
