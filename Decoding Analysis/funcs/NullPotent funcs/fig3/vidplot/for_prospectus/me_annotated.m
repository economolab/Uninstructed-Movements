close all

data = me.data(:,params.trialid{1});
move = me.move(:,params.trialid{1});

ixs = [180 240];
[~,sortix] = sort(mean(data(ixs(1):ixs(2),:),1),'descend');

n = 25; 
trix = randsample(sortix(1:200),n,false);

data = data(:,trix);
move = move(:,trix);


%%
close all
tempdat = data;
tempdat = mySmooth(tempdat,11);
figure;
ax = gca;
hold on;
dy = 50;
for i = 1:n
    temp = tempdat(:,i);
%     temp(move(:,i)) = nan;
    plot(obj.time, temp + (i)*dy,'k','LineWidth',1.5)
    temp = tempdat(:,i);
    temp(~move(:,i)) = nan;
    plot(obj.time, temp + (i)*dy,'r','LineWidth',1.5)
end
xs = [-1 -0.02 -0.02 -1];
ys = [ax.YLim(1) ax.YLim(1) ax.YLim(2), ax.YLim(2)];
fill(xs,ys,'k','EdgeColor','none','FaceAlpha',0.2)

xs = [0.02 1 1 0.01];
ys = [ax.YLim(1) ax.YLim(1) ax.YLim(2), ax.YLim(2)];
fill(xs,ys,'k','EdgeColor','none','FaceAlpha',0.5)

xlim([-2.4 2.5])
ylim([ax.YLim(1) 1350])














