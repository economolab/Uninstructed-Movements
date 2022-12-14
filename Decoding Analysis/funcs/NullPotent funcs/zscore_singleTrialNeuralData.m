function trialdat_zscored = zscore_singleTrialNeuralData(dat)
% zscore obj.trialdat, which is single trial binned neural data
% (time,neurons,trials)
% The zscoring takes place on the matrix of (time*trials,neurons) for all
% neurons
% trialdat_zscored is of size(time,trials,neurons)



temp = permute(dat,[1 3 2]);
dims = size(temp); % (time,trials,neurons)

temp = reshape(temp,dims(1)*dims(2),dims(3));

temp2 = zscore(temp);
trialdat_zscored = reshape(temp2,dims(1),dims(2),dims(3));




end