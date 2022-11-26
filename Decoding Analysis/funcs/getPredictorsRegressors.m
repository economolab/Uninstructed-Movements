function [Y,X] = getPredictorsRegressors(trials,regr,kin,rez)
% OUTPUT: What you are trying to predict--CDlate for all of the train trials
% Y = time x trials
Y = [];
for t = 1:length(trials.all)
    currtrix = trials.all(t);
    Y = [Y,regr.singleProj(:,currtrix)];
end

% INPUT: What you are using to predict--kinematic measures for the given body part for all of the trian trials
% X = time x trials x num measures for body part
X = kin.dat(:,trials.all,rez.featix);
% fill missing values
for featix = 1:size(X,3)
    X(:,:,featix) = fillmissing(X(:,:,featix),"constant",0);
end
end