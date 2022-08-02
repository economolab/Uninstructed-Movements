function acc = train_test_choiceDecoder(acc,X_train,X_test,y_train,y_test,bootiter,dt,sessix)
ix = 1;
for i = 1:dt:size(X_train,1) % each timepoint
    x_train = squeeze(X_train(i,:,:)); % (trials,feats) for timepoint i

    if size(x_train,2) > 1 && size(x_train,1) ~= 1
        mdl = fitcecoc(x_train,y_train); % uses binary svm classifier by default
    else
        mdl = fitcecoc(x_train',y_train);
    end

    x_test = squeeze(X_test(i,:,:)); % (trials,feats) for timepoint i
    if size(x_test,2) > 1 && size(x_test,1) ~= 1
        pred = predict(mdl,x_test);
    else 
        pred = predict(mdl,x_test');
    end
    acc(ix,bootiter,sessix) = sum(pred == y_test) / numel(y_test);
    ix = ix + 1;
end
end