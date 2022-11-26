dat = me(1);
obj_ = obj(1);

k = 20; % trials to sample
nTrials = size(dat.data,2);

trix = randsample(nTrials,k,'false');

%%

close all

sm = 7;
dy = 60;

f = figure; 
hold on
for i = 1:numel(trix)
    t = trix(i);
    temp = mySmooth(dat.data(:,t),sm);
    move = dat.move(:,t);
    temp2 = temp;
    temp2(~move) = nan;
    plot(obj_.time,temp + (i*dy),'Color',[150 150 150]./255,'LineWidth',1.5)
    plot(obj_.time,temp2 + (i*dy),'r','LineWidth',1.5)
end

align = mode(obj_.bp.ev.goCue);
sample = mode(obj_.bp.ev.sample) - align;
delay = mode(obj_.bp.ev.delay) - align;
xline(sample,'k--','LineWidth',1)
xline(delay,'k--','LineWidth',1)
xline(0,'k--','LineWidth',1)

xlim([-2.4, 2])
xlabel('Time (s) from go cue')
ylabel('Trials')
title('Motion Energy | JEB15 - 2022-07-26')
ax = gca;
ax.FontSize = 14;
ax.YTick = [];


%%
close all

data = me.data;

z = data - (me.moveThresh+5);
z(z<=0) = -max(max(z));

f = figure; hold on;
f.Position = [-1311        -179         798        1003];
imagesc(obj.time,1:size(z,2),z')

xline(sample,'k--','LineWidth',3)
xline(delay,'k--','LineWidth',3)
xline(0,'k--','LineWidth',3)

hold off

ylim([1,size(z,2)])

a = colorbar;
a.Label.String = '(Motion Energy - Threshold) > 0';
xlabel('Time (s) from go cue')
ylabel('Trials')
ax = gca;
ax.FontSize = 30;

colormap cool
% 
% pth = '/Users/Munib/Documents/Economo-Lab/code/uninstructedMovements/fig4/figs/motionEnergy/';
% fn = ['motionEnergyMap_JEB7_2021-04-29_threshincreasedby5'];
% mysavefig(f,pth,fn);