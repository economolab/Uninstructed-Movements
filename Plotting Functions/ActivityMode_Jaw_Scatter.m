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
        jawToUse = jawVel(cnt+1:nTrials);
    else
        jawToUse = jawVel(cnt:nTrials); 
    end
    cnt = cnt+1+nTrials;

    if strcmp(params.earlytrials,'motionEnergy')
        temp = find(obj.earlyMoveix);
        %earlycond = ismember(trialid,temp)
    else
        temp = obj.earlyMoveix;
    end
    earlyTrialid = ismember(trialid,temp);
    regTrialid = ~ismember(trialid,temp);
    color = clrs{j};                     
    scatter(jawToUse(regTrialid),avgChoice(regTrialid),'MarkerFaceColor',color)     % Make a scatter plot comparing these values 
    hold on;
    scatter(jawToUse(earlyTrialid),avgChoice(earlyTrialid),'MarkerFaceColor',color,'Marker','diamond')                                                                        % for the desired trials, colored according to condition
end
ax = gca;
ax.FontSize = 12;
end  %ActivityMode_Jaw_Scatter
