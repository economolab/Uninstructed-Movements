function plotFeatVelocity(traj,vel,trix)

time = traj(trix).time;

[datmat,labels] = featStruct2Mat(vel(trix));

f = figure(101);
imagesc('XData',time,'YData',1:numel(labels),'CData',datmat')
clim1 = min(min(datmat))/5;
clim2 = max(max(datmat))/5;
caxis([clim1 clim2])
xlim([time(1) time(end)])
ylim([1 numel(labels)])
xlabel('Time (s) from movement onset')
ylabel('Features')
title(['Trial ' num2str(trix)])
ax = f.CurrentAxes;
ax.YTickLabel = labels;
ax.FontSize = 20;
colormap hot
cbar = colorbar;
cbar.Label.String = 'Velocity';
set(f,'Position',[-1915         224        1907         565])

end












