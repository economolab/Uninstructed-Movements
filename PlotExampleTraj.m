figure();
temp=obj.traj{1};
trial = 5;
prevmax = 0;
allLicks = [obj.bp.ev.lickR{trial},obj.bp.ev.lickL{trial}]; allLicks = allLicks+0.5;
allLicks = sort(allLicks,'ascend');

for i=1:numel(temp(1).featNames)
    toplot = temp(trial).ts(:,2,i);
    plot(temp(trial).frameTimes, prevmax/2+toplot,'LineWidth',1.5)
    prevmax = prevmax+max(toplot);
    hold on;
end
ylicks = (prevmax/2)+10*ones(1,length(allLicks));
scatter(allLicks,ylicks,'black','filled')
xline(obj.bp.ev.goCue(trial)+0.5,'LineStyle','--')
xline(obj.bp.ev.delay(trial)+0.5,'LineStyle','-.')
xline(obj.bp.ev.sample(trial)+0.5,'LineStyle','-.')
legend(temp(1).featNames,'Location','best')
trial = num2str(trial);
trialtitle = strcat('Example traj for',{' '},anm,{' '},date,';',{' '},'Trial',{' '},trial);  % Name/title for session
title(trialtitle)
