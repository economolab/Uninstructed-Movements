function R2 = plotNP_PredictedCDLate_Scatter(decoded, sessix, meta, string,start,stop)
tempR = mean(decoded(sessix).trueVals.Rhit((start+1):(stop+1),:),1,'omitnan');        % For each trial, get the average CDlate during the delay period
    tempL = mean(decoded(sessix).trueVals.Lhit((start+1):(stop+1),:),1,'omitnan');
    truedat = [tempR, tempL];

    tempR = mean(decoded(sessix).modelpred.Rhit(start:stop,:),1,'omitnan');        % For each trial, get the average predicted CDlate during the delay period
    tempL = mean(decoded(sessix).modelpred.Lhit(start:stop,:),1,'omitnan');
    modeldat = [tempR, tempL];

    R2 = corrcoef(truedat,modeldat);
    R2 = R2(2);
    sR2 = num2str(R2);
    coeff = polyfit(truedat,modeldat,1);                   % Find the line of best fit
    sesstitle = [meta(sessix).anm, ' ',meta(sessix).date];
    if strcmp(string,'null')
        sub = 1;
    elseif strcmp(string,'potent')
        sub = 2;
    end
    subplot(1,2,sub)
    scatter(truedat,modeldat,15,'filled','black')
    hline = refline(coeff);
    hline.LineStyle = '--';
    hline.Color = 'k';
    xlabel('True data')
    ylabel('Model prediction')
    legend('data',['R^2 = ' sR2],'Location','best')
    title(string)
    
    sgtitle(sesstitle)
end