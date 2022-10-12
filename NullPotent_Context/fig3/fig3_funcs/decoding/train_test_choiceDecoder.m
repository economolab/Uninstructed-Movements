function acc = train_test_choiceDecoder(acc,X_train,X_test,y_train,y_test,bootiter,binSize)

ix = 1;
for i = 1:binSize:(size(X_train,1)-binSize) % each timepoint

    ixs = i:(i+binSize-1);

    % train
    x_train = X_train(ixs,:);

    mdl = fitcecoc(x_train',y_train);

    % test
    x_test = X_test(ixs,:);

    pred = predict(mdl,x_test');
    
    acc(ix,bootiter) = sum(pred == y_test) / numel(y_test);
    ix = ix + 1;
end
end






