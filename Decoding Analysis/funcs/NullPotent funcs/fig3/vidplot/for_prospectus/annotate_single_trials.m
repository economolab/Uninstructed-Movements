

trialdat = obj.trialdat(:,:,params.trialid{1});

clu = randsample(size(trialdat,2),3,false);
tempdat = trialdat(:,clu,trix);


%%
close all
dy = 1.25;
for i = 1:numel(clu)
    figure;
    ax = gca;
    hold on;
    for j = 1:numel(trix)
        temp = mySmooth(normalize(tempdat(:,i,j),'range',[0 1]),11);
        plot(obj.time, temp + (j)*dy,'k','LineWidth',1.5)
        temp = mySmooth(normalize(tempdat(:,i,j),'range',[0 1]),11);
        temp(~move(:,j)) = nan;
        plot(obj.time, temp + (j)*dy,'r','LineWidth',1.5)
    end
    xlim([-2.4 2.5])
    ylim([ax.YLim(1) 33])
end

