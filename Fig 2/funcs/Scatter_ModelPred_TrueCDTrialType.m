function [R2] = Scatter_ModelPred_TrueCDTrialType(trueVals, modelpred, sessix, start, stop,meta,invert)
if strcmp(invert,'invert')
    trueR = mean(trueVals.Rhit{sessix}(start:stop,:),1,'omitnan'); trueR = -1*trueR;        % For each trial, get the average CDlate during the delay period
    trueL = mean(trueVals.Lhit{sessix}(start:stop,:),1,'omitnan'); trueL = -1*trueL;
    truedat = [trueR, trueL];

    modR = mean(modelpred.Rhit{sessix}(start:stop,:),1,'omitnan'); modR = -1*modR;       % For each trial, get the average predicted CDlate during the delay period
    modL = mean(modelpred.Lhit{sessix}(start:stop,:),1,'omitnan'); modL = -1*modL; 
    modeldat = [modR, modL];
else
    trueR = mean(trueVals.Rhit{sessix}(start:stop,:),1,'omitnan');        % For each trial, get the average CDlate during the delay period
    trueL = mean(trueVals.Lhit{sessix}(start:stop,:),1,'omitnan');
    truedat = [trueR, trueL];

    modR = mean(modelpred.Rhit{sessix}(start:stop,:),1,'omitnan');        % For each trial, get the average predicted CDlate during the delay period
    modL = mean(modelpred.Lhit{sessix}(start:stop,:),1,'omitnan');
    modeldat = [modR, modL];
end

R2 = corrcoef(truedat,modeldat);
R2 = R2(2);
coeff = polyfit(modeldat,truedat,1);                   % Find the line of best fit
sesstitle = strcat(meta(sessix).anm, {' '},meta(sessix).date);
%scatter(modeldat,truedat,20,'filled','black')
scatter(modL,trueL,20,'filled','red'); hold on
scatter(modR,trueR,20,'filled','blue')
hline = refline(coeff);
hline.LineStyle = '--';
hline.Color = 'k';
xlabel('Model prediction')
ylabel('True data')
legend('data',['R^2 = ' num2str(R2)],'Location','Best')
title(sesstitle)
