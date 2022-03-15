% Plots a scatter plot of jaw velocity vs. choice mode on individual trials
% Each data point = an average value of jaw vel or activity mode during a
% specific time epoch

% INPUTS: jawVel = [1 x nTrials] array of average jaw velocities
% avgChoice = [1 x nTrials] array of average choice mode values
% conditions = range of numbers to indicate which conditions you want to
% make this plot for (i.e. 1:2 will plot R and L 2AFC hits)
% clrs = cell array of colors to use

function [coeff,R] = ActivityMode_Jaw_Scatter(jawVel,avgChoice,conditions,met,clrs,obj,params)
jv = [];                                        
ch = [];
for j = conditions                              % For all conditions                 
    trialid = met.trialid{j};                   % Get the trials that correspond to those conditions
    if strcmp(params.earlytrials,'motionEnergy')
        temp = find(obj.earlyMoveix);
    else
        temp = obj.earlyMoveix;
    end
    earlyTrialid = trialid(find(ismember(trialid,temp)));
    regTrialid = trialid(find(~ismember(trialid,temp)));
    color = clrs{j};                    
    jv = [jv, jawVel(trialid)];                 % Store only the jaw velocities from the desired trials 
    ch = [ch, avgChoice(trialid)];              % Store only the choice mode values from the desired trials 
    scatter(jawVel(regTrialid),avgChoice(regTrialid),'MarkerFaceColor',color)     % Make a scatter plot comparing these values 
    hold on;
    scatter(jawVel(earlyTrialid),avgChoice(earlyTrialid),'MarkerFaceColor',color,'Marker','diamond')                                                                        % for the desired trials, colored according to condition
end
ax = gca;
ax.FontSize = 12;
nanix = find(isnan(jv));                    % Find indices where the jaw vel is a NaN
jv = jv(~isnan(jv));                        % Get rid of the NaN values
ch(nanix) = [];                             % Indices that were a NaN for jaw vel, get rid of those indices in the choice mode as well
R = corr2(jv,ch);                           % Calculate the correlation coefficient between these two variables
R = num2str(R);
coeff = polyfit(jv,ch,1);                   % Find the line of best fit
end  %ActivityMode_Jaw_Scatter
