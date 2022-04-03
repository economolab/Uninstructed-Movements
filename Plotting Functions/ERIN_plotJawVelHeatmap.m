function ERIN_plotJawVelHeatmap(alljaw)
toplot = [];
taxis = 1:1401;
for i = 1:numel(alljaw)
    toplot = [toplot,alljaw{i}];
end
l1 = size(alljaw{1},2);
l2 = l1+size(alljaw{2},2);
l3 = l2+size(alljaw{3},2);
imagesc(toplot')
colorbar
caxis([0 1.5])
hold on;
line([taxis(1),taxis(end)],[l1,l1],'Color','white','LineStyle','--')
line([taxis(1),taxis(end)],[l2,l2],'Color','white','LineStyle','--')
line([taxis(1),taxis(end)],[l3,l3],'Color','white','LineStyle','--')
end