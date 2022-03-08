function plotKMeansEthogram(trix,time,datmat,labels,idx)


f = figure(101);
subplot(211)
imagesc('XData',time,'YData',1:numel(labels),'CData',datmat')
clim1 = min(min(datmat))/5;
clim2 = max(max(datmat))/5;
caxis([clim1 clim2])
xlim([time(1) time(end)])
ylim([1 numel(labels)])
ylabel('Features')
title(['Trial ' num2str(trix)])
ax = f.CurrentAxes;
ax.YTickLabel = labels;
ax.FontSize = 20;
colormap hot
cbar = colorbar;
cbar.Location = 'southoutside';
cbar.Label.String = 'Velocity';

subplot(212)
y = 1;
imagesc('XData',time,'YData',y,'CData',idx')
xlim([time(1) time(end)])
ylim([y-0.5 y+0.5])
xlabel('Time (s) from movement onset')
ylabel('Ethogram')
ax = f.CurrentAxes;
ax.YTick = [];
ax.FontSize = 20;
ax.Colormap = parula;
set(f,'Position',[-1915         224        1907         565])


end












