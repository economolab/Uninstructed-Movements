function [R2] = Scatter_ModelPred_TrueCDTrialType(trueVals, modelpred, sessix, start, stop,meta)
tempR = mean(trueVals.Rhit{sessix}((start+1):(stop+1),:),1,'omitnan');        % For each trial, get the average CDlate during the delay period
tempL = mean(trueVals.Lhit{sessix}((start+1):(stop+1),:),1,'omitnan');
truedat = [tempR, tempL];

tempR = mean(modelpred.Rhit{sessix}(start:stop,:),1,'omitnan');        % For each trial, get the average predicted CDlate during the delay period
tempL = mean(modelpred.Lhit{sessix}(start:stop,:),1,'omitnan');
modeldat = [tempR, tempL];

R2 = corrcoef(truedat,modeldat);
R2 = R2(2);
coeff = polyfit(modeldat,truedat,1);                   % Find the line of best fit
sesstitle = strcat(meta(sessix).anm, {' '},meta(sessix).date);
scatter(modeldat,truedat,20,'filled','black')
hline = refline(coeff);
hline.LineStyle = '--';
hline.Color = 'k';
xlabel('Model prediction')
ylabel('True data')
legend('data',['R^2 = ' num2str(R2)],'Location','Best')
title(sesstitle)
