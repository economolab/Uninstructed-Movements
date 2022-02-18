% Plots a scatter plot of jaw velocity vs. choice mode on individual trials
% Each data point = an average value of jaw vel or activity mode during a
% specific time epoch

% INPUTS: jawVel = [1 x nTrials] array of average jaw velocities
% avgChoice = [1 x nTrials] array of average choice mode values
% conditions = range of numbers to indicate which conditions you want to
% make this plot for (i.e. 1:2 will plot R and L 2AFC hits)
% clrs = cell array of colors to use

function ActivityMode_Jaw_Scatter(jawVel,avgChoice,conditions,met,clrs)
jv = [];                                        
ch = [];
for j = conditions                              % For all conditions                 
    trialid = met.trialid{j};                   % Get the trials that correspond to those conditions
    color = clrs{j};                    
    jv = [jv, jawVel(trialid)];                 % Store only the jaw velocities from the desired trials 
    ch = [ch, avgChoice(trialid)];              % Store only the choice mode values from the desired trials 
    scatter(jawVel(trialid),avgChoice(trialid),'MarkerFaceColor',color)     % Make a scatter plot comparing these values 
                                                                            % for the desired trials, colored according to condition
    hold on;
end
ax = gca;
ax.FontSize = 12;
nanix = find(isnan(jv));                    % Find indices where the jaw vel is a NaN
jv = jv(~isnan(jv));                        % Get rid of the NaN values
ch(nanix) = [];                             % Indices that were a NaN for jaw vel, get rid of those indices in the choice mode as well
R = corr2(jv,ch);                           % Calculate the correlation coefficient between these two variables
R = num2str(R);
coeff = polyfit(jv,ch,1);                   % Find the line of best fit
hline = refline(coeff);
hline.LineStyle = '--';
hline.Color = 'k';
str = strcat('R^2 =',R);
title(str,'fontsize',13)
xlabel('Avg Jaw Velocity','fontsize',14)
ylabel('Avg choice mode','fontsize',14)
end  %ActivityMode_Jaw_Scatter
