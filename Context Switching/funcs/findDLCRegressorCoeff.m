function beta = findDLCRegressorCoeff(in,rez)

% bootstrap, and sessions happen outside this function
% This function is for a single session, single iteration
X_train = in.train.X;
Y_train = in.train.y;
X_test = in.test.X;
Y_test = in.test.y;

nRegressors = size(X_train,3);

TimePtsInTrain = size(X_train,1);       % Number of time points in the training data
nTimePts = TimePtsInTrain - rez.dt;     % Have to subtract the number of samples in the time bin (cannot use the last time point because will not be able to create the bin)
tRange = 1:rez.dt:nTimePts;             % Will loop over all of these time-points 

ix = 1;
for i = tRange                          % For each time-point...

    ixs = i:(i+rez.dt-1);               % Get the time indices that are in the current bin for this time-point

    % Get regressor and predictor data for the current time bin
    % 'x_train' ends up being (trials x features) 
    x_train = X_train(ixs,:,:);       % Get the regressor training data for the current time-bin 
    x_train = squeeze(mean(x_train,1,'omitnan'));   % Take the average value of the regressor across time 
    x_train = fillmissing(x_train,'constant',0);    % Make NaNs into zeros
    outlierx = find(abs(x_train)>100);              % If any of the regressors (movement values) have a very large value (greater than 100), label it as an outlier
    if ~isempty(outlierx)                   
        x_train(outlierx) = 0;                      % Set these outlier values to 0
    end
    if size(x_train,1) == 1
        x_train = x_train';
    end

    % 'y_train' ends up being (trials x 1)
    y_train = Y_train(ixs,:,:);                     % Get the predictor training data for the current time-bin
    y_train = squeeze(mean(y_train,1,'omitnan'));
    y_train = fillmissing(y_train,'constant',0);    % Make NaNs into zeros
    if size(y_train,1) == 1
        y_train = y_train';
    end

    % Fit the linear regression model
    %cv_mdl = fitrlinear(x_train,y_train);
    cv_mdl = fitrlinear(x_train,y_train,'KFold',rez.nFolds);    % Use cross-validation to find model

    tempbet = NaN(nRegressors,rez.nFolds);
    for k = 1:rez.nFolds                                        % For each iteration of cross-validation, take the Beta vals ('weight' for each movement feature)
        tempbet(:,k) = cv_mdl.Trained{k}.Beta;
    end
    beta(:,ix) = mean(tempbet,2,'omitnan');                     % Average them across CV iterations, and store this for the time-point

    ix = ix+1;
end

end