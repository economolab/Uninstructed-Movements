figure();
temp=obj.traj{1};
trial = 5;
prevmax = 0;
allLicks = [obj.bp.ev.lickR{trial},obj.bp.ev.lickL{trial}]; allLicks = allLicks+0.5;
allLicks = sort(allLicks,'ascend');

[ang,len] = findTongueAngle(obj,trial);

for i=[2,4,6,8]
    if i~=8
        toplot = temp(trial).ts(:,2,i);
        plot(temp(trial).frameTimes, (prevmax/2)+toplot,'LineWidth',1.5,'LineStyle','-')
        prevmax = prevmax+max(toplot);
        hold on;
    elseif i==8
        plot(temp(trial).frameTimes,300+(mySmooth(me.data{trial},200)),'LineWidth',1.5,'LineStyle','-')
    end
end

ylicks = (prevmax-120)*ones(1,length(allLicks));
scatter(allLicks,ylicks,'black','filled')
xline(obj.bp.ev.goCue(trial)+0.5,'LineStyle','--')
xline(obj.bp.ev.delay(trial)+0.5,'LineStyle','-.')
xline(obj.bp.ev.sample(trial)+0.5,'LineStyle','-.')
xlim([0.5 6])
legend('Tongue','Jaw','Nose','Motion Energy','Location','best')
trial = num2str(trial);
anm = obj.meta.anm;
date = obj.meta.day;
trialtitle = strcat('Example traj for',{' '},anm,{' '},date,';',{' '},'Trial',{' '},trial);  % Name/title for session
title(trialtitle)

%%
function [ang,len] = findTongueAngle(obj,trial)
dat = medfilt1(obj.traj{2}(trial).ts, 3, [], 1);
    
%Find the median x and y jaw position for the trial i
jx = nanmedian(dat(:, 1, 8));
jy = nanmedian(dat(:, 2, 8));

%Find the x and y tongue tip position for the all time points in trial i
tx = (dat(:, 1, 1)+dat(:, 1, 3))./2;    %Average between x position of top tongue and bottom tongue
ty = (dat(:, 2, 1)+dat(:, 2, 3))./2;    %Average between y position of top tongue and bottom tongue

dx = tx-jx;                             %Distance in x coordinates between tongue tip and jaw
dy = ty-jy;                             %Distance in y coordinates between tongue tip and jaw
len = sqrt(dx.^2 + dy.^2);           %Length of tongue for all points in trial i
ang = atan(dy./dx);                  %Angle of tongue for all points in trial i

ang(dx<0 & dy>0) = ang(dx<0 & dy>0) + pi;     %Correction for the quadrant that the angle lies in
ang(dx<0 & dy<0) = ang(dx<0 & dy<0) - pi;
ang = mySmooth(ang,4);
end