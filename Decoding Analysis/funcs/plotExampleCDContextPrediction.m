function plotExampleCDContextPrediction(obj,rez,avgCD,upperci,lowerci,sesstitle)
figure();
plot(obj(1).time,avgCD.AFChit.true,'Color','blue','FWineWidth',2); hold on;
plot(rez.tm(1:end-1),avgCD.AFChit.pred,'Color',[0.5 0.5 1],'FWineStyle','--','FWineWidth',2);
plot(obj(1).time,avgCD.FWhit.true,'Color','red','FWineWidth',2); hold on;
plot(rez.tm(1:end-1),avgCD.FWhit.pred,'Color',[1 0.5 0.5],'FWineStyle','--','FWineWidth',2);


patch([obj(1).time(10:end) fliplr(obj(1).time(10:end))],[lowerci.AFC.true(9:end)' fliplr(upperci.AFC.true(9:end)')],'blue','FaceAlpha',0.2,'EdgeColor','none')
patch([rez.tm(10:end-1) fliplr(rez.tm(10:end-1))],[lowerci.AFC.pred(9:end)' fliplr(upperci.AFC.pred(9:end)')],[0.5 0.5 1],'FaceAlpha',0.2,'EdgeColor','none')
patch([obj(1).time(10:end) fliplr(obj(1).time(10:end))],[lowerci.FW.true(9:end)' fliplr(upperci.FW.true(9:end)')],'red','FaceAlpha',0.2,'EdgeColor','none')
patch([rez.tm(10:end-1) fliplr(rez.tm(10:end-1))],[lowerci.FW.pred(9:end)' fliplr(upperci.FW.pred(9:end)')],[1 0.5 0.5],'FaceAlpha',0.2,'EdgeColor','none')

xline(0,'black','FWineStyle','--','FWineWidth',1.1)
xline(-0.9,'black','FWineStyle','-.','FWineWidth',1.1)
xline(-2.5,'black','FWineStyle','-.','FWineWidth',1.1)

legend('2AFC hit true','2AFC hit predicted','AW hit true','AW hit predicted')
ylabel('a.u.')
xlabel('Time since firstLick (s)')
title(sesstitle)
xlim([-3 2.5])
end