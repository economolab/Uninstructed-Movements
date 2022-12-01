function plotExampleCDContextPrediction(obj,rez,avgCD,upperci,lowerci,sesstitle)
figure();
plot(obj(1).time,avgCD.AFChit.true,'Color','black','LineWidth',2); hold on;
plot(rez.tm(1:end-1),avgCD.AFChit.pred,'Color','black','LineStyle','--','LineWidth',2);
plot(obj(1).time,avgCD.FWhit.true,'Color','magenta','LineWidth',2); hold on;
plot(rez.tm(1:end-1),avgCD.FWhit.pred,'Color','magenta','LineStyle','--','LineWidth',2);


patch([obj(1).time(10:end) fliplr(obj(1).time(10:end))],[lowerci.AFC.true(9:end)' fliplr(upperci.AFC.true(9:end)')],'black','FaceAlpha',0.2,'EdgeColor','none')
patch([rez.tm(10:end-1) fliplr(rez.tm(10:end-1))],[lowerci.AFC.pred(9:end)' fliplr(upperci.AFC.pred(9:end)')],'black','FaceAlpha',0.2,'EdgeColor','none')
patch([obj(1).time(10:end) fliplr(obj(1).time(10:end))],[lowerci.FW.true(9:end)' fliplr(upperci.FW.true(9:end)')],'magenta','FaceAlpha',0.2,'EdgeColor','none')
patch([rez.tm(10:end-1) fliplr(rez.tm(10:end-1))],[lowerci.FW.pred(9:end)' fliplr(upperci.FW.pred(9:end)')],'magenta','FaceAlpha',0.2,'EdgeColor','none')

xline(0,'black','LineStyle','--','LineWidth',1.1)
xline(-0.9,'black','LineStyle','-.','LineWidth',1.1)
samp = median(obj(1).bp.ev.sample)-median(obj(1).bp.ev.goCue);
xline(samp,'black','LineStyle','-.','LineWidth',1.1)

legend('2AFC hit true','2AFC hit predicted','AW hit true','AW hit predicted','Location','best')
ylabel('a.u.')
xlabel('Time since goCue (s)')
title(sesstitle)
xlim([-3 2.5])
end