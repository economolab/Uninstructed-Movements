function acc = DLC_ChoiceDecoder(in,rez,trials)

% bootstrap, and sessions happen outside this function

X_train = in.train.X;
y_train = in.train.y;
X_test = in.test.X;
y_test = in.test.y;


ix = 1;
for i = 1:rez.dt:(size(X_train,1)-rez.dt) % each timepoint

    ixs = i:(i+rez.dt-1);

    % train
    x_train = squeeze(nanmean(X_train(ixs,:,:),1));
    x_train = fillmissing(x_train,'constant',0);
    if size(x_train,1) == 1
        x_train = x_train';
    end

%     mdl = fitcecoc(x_train,y_train);
    mdl = fitclinear(x_train',y_train','Learner','logistic','KFold',rez.nFolds,'ObservationsIn','columns');
    % mdl = fitcsvm(x_train,y_train,'Standardize',false);
    % cv_mdl = crossval(mdl,'KFold',rez.nFolds);
    
    pred = kfoldPredict(mdl);

    acc(ix) = sum(pred == y_train) / numel(y_train);
    ix = ix + 1;
end


end