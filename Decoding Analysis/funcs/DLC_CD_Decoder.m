function pred = DLC_CD_Decoder(in,rez)

% bootstrap, and sessions happen outside this function
% This function is for a single session, single iteration
X_train = in.train.X;
Y_train = in.train.y;
X_test = in.test.X;
Y_test = in.test.y;

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
    if size(x_train,1) == 1
        x_train = x_train';
    end

    % 'y_train' ends up being (trials x 1)
    y_train = Y_train(ixs,:,:);
    y_train = squeeze(mean(y_train,1,'omitnan'));
    y_train = fillmissing(y_train,'constant',0);    % Make NaNs into zeros
    if size(y_train,1) == 1
        y_train = y_train';
    end

    % Fit the linear regression model
    mdl = fitrlinear(x_train,y_train);              
    % cv_mdl = crossval(mdl,'KFold',rez.nFolds);


    % Test the model on test data 
    % Get the average test data for the given time-bin
    x_test = squeeze(mean(X_test(ixs,:,:),1,'omitnan'));
    x_test = fillmissing(x_test,'nearest');
    if size(x_test,1) == 1
        x_test = x_test';
    end
    
    % Will predict CDlate for the current time-bin across all trials
    pred(:,ix) = predict(mdl,x_test);       
%     pred = kfoldPredict(cv_mdl);
%     pred = predict(cv_mdl,x_test);

%     acc(ix) = sum(pred == y_test) / numel(y_test);
%    acc(ix) = sum(pred == y_train) / numel(y_train);
    ix = ix + 1;
end

end