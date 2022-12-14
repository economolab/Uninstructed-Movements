function acc = choiceDecoder(dat,numT,k,params,train,binSize,cond2use)

acc = zeros(numT,k); % (time,bootstraps)


for bootix = 1:k
    disp(['--Iteration ' num2str(bootix) '/' num2str(k)])

    y_train = [];
    y_test = [];
    X_train = [];
    X_test = [];
    for condix = 1:numel(cond2use) % looping throgh conditions to build data
        % get trials for current session and condition
        trials2use = params.trialid{cond2use(condix)};
        nTrials = numel(trials2use);

        % number of train and test trials
        nTrain = round(train*nTrials); % per condition
        nTest  = nTrials - nTrain;     % per condition

        % get training and testing trials
        train_trials = randsample(trials2use,nTrain);
        remain = trials2use(~ismember(trials2use,train_trials));
        test_trials = randsample(remain,nTest);

        y_train_cond = condix*ones(size(train_trials));
        y_test_cond = condix*ones(size(test_trials));

        X_train_cond = dat(:,train_trials);
        X_test_cond = dat(:,test_trials);

        y_train = cat(1,y_train,y_train_cond);
        y_test = cat(1,y_test,y_test_cond);

        X_train = cat(2,X_train,X_train_cond);
        X_test = cat(2,X_test,X_test_cond);

    end

    y_train(y_train==2) = -1;
    y_test(y_test==2) = -1;

    % choice decoding

    acc = train_test_choiceDecoder(acc,X_train,X_test,...
        y_train,y_test,bootix,binSize);

end



end