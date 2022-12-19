function R2 = correlateTrue_PredictedCDContext(trueVals, modelpred,meta)
% This is across sessions
R2 = NaN(length(meta),1);

for sessix = 1:length(meta)                         % For each session...
tempR1 = [];
for t = 1:size(trueVals.AFChit{sessix},2)                 % For every 2AFc hit trial...
    trueR = trueVals.AFChit{sessix}(1:end-1,t);   % Take the true CDlate value
    predR = modelpred.AFChit{sessix}(:,t);        % Take the predicted CDlate value
    predR = fillmissing(predR,"nearest");       % Get rid of NaNs
    corr = corrcoef(trueR,predR);               % Find the correlation between the true and predicted CDlate
    tempR1 = [tempR1;corr(1,2)];                % Store this R2 value for each trial
end
tempR2 = [];
for t = 1:size(trueVals.FWhit{sessix},2)         % Do the same for AW hit trials
    trueL = trueVals.FWhit{sessix}(1:end-1,t);
    predL = modelpred.FWhit{sessix}(:,t);
    predL = fillmissing(predL,"nearest");
    corr = corrcoef(trueL,predL);
    tempR2 = [tempR2;corr(1,2)];
end

R2(sessix) = mean([tempR1;tempR2],'omitnan');   % Take the average R2 value across all trials for this session
end
end