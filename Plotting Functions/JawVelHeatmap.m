function JawVelHeatmap(conditions, jaw, taxis, met)
toplot = [];
for c = 1:numel(conditions)
     toplot = [toplot; met.trialid{c}];  % Trial ID of right, hit, non-early 2AFC trials
end

trix = sort(toplot);                    % Trial ID of all hit, non-early 2AFC trials
numtrix = length(trix);              % Num of trials to be plotted

imagesc(taxis(10:end), 1:numtrix,jaw(10:end,trix)')      % Make heatmap of jaw velocity for specified trials
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