% Plots a heatmap of jaw velocities for all trials in a session

function JawVelHeatmap(conditions, jaw, taxis)
temp = [];
for c = 1:numel(conditions)
    temp = [temp,jaw{c}];
end
numtrix = size(temp,2);

imagesc(taxis(10:end), 1:numtrix,temp(10:end,:)')      % Make heatmap of jaw velocity for specified trials
ax = gca;
ax.FontSize = 12;
c = colorbar;
lim = [0 7];
caxis(lim)
ylabel(c,'Velocity')
xlabel('Time since go-cue (s)','FontSize',13)
ylabel('Trials','FontSize',13)
title('Heatmap of jaw velocity','FontSize',14)
end   % JawVelHeatmap