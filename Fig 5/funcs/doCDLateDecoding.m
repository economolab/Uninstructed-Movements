function [trueVals, modelpred] = doCDLateDecoding(params,sessix,cond2use,hitcond,misscond,space,kin,rez,string)

% Getting the proper number of trials
[trials,trials_hit] = EqualizeTrialNums(params(sessix),cond2use,hitcond,misscond);

% Organize predictor and regressors from the proper trials and
% kinematic features
[Y,X] = getPredictorsRegressors_NullPotent(trials,space(sessix),kin(sessix),rez);

% Train/Test Split
[trials,in] = TrainTestSplit(trials,Y,X,rez);

% Decoding
pred = DLC_CD_Decoder(in,rez);
pred = mySmooth((pred'),31);

% Divide true data and model predictions into R and L hit trials
if rez.train==1                                                     % If all data is used to train the model (Cross-validated, the true data is just all of it i.e. the training data)
    trials.RHit.TrainIX = ismember(trials.train,trials_hit{1});
    trials.RHit.Train = trials.train(trials.RHit.TrainIX);
    trials.LHit.TrainIX = ismember(trials.train,trials_hit{2});
    trials.LHit.Train = trials.train(trials.LHit.TrainIX);
else                                                                % If model is not cross-validated and just doing a train/test split
    trials.RHit.TestIX = ismember(trials.test,trials_hit{1});       % Logical array for all of the test trials indicating whether they were a right hit or not
    trials.RHit.Test = trials.test(trials.RHit.TestIX);             % Get the trial numbers that are a R hit and were used in the test
    trials.LHit.TestIX = ismember(trials.test,trials_hit{2});       % Logical array for all of the test trials indicating whether they were a left hit or not
    trials.LHit.Test = trials.test(trials.LHit.TestIX);             % Get the trial numbers that are a L hit and were used in the test
end

% Label test data as true data and model prediction
if rez.train==1                                                     % If you are using the whole data set as the training set (if cross-validating), then the test data is just the train data
    trueVals.Rhit = in.train.y(:,trials.RHit.TrainIX);
    trueVals.Lhit= in.train.y(:,trials.LHit.TrainIX);
    modelpred.Rhit = pred(:,trials.RHit.TrainIX);
    modelpred.Lhit = pred(:,trials.LHit.TrainIX);
else                                                                % If you are using a train/test split, then the assign the test data as the data not used for training
    trueVals.Rhit = in.test.y(:,trials.RHit.TestIX);
    trueVals.Lhit = in.test.y(:,trials.LHit.TestIX);
    modelpred.Rhit = pred(:,trials.RHit.TestIX);
    modelpred.Lhit = pred(:,trials.LHit.TestIX);
end
end