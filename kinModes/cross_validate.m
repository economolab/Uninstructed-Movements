function best_lambda = cross_validate(Y,X,lambdas)
n = size(X,1); % num observations (time)
K = 5; % fraction of data held out for testing
c = cvpartition(n,'Kfold',K);
% for each alpha
%   for each k in k-fold
%       keep fold k as hold-out data
%       use remaining folds and current alpha to estimate W
%       predict held-out data: M_test,k = X_test,k * W
%       compute MSE: |M - M_test,k|^2
%   end for k
%   average MSE over the k folds: 1/K * sum(MSE)
% end for alpha
% choose optimal value: alpha_opt = argmin_p (MSE_p)

mse_i = zeros(1,length(lambdas));
for i = 1:length(lambdas)
    mse_k = zeros(1,K);
    for k = 1:K
        idxTrain = training(c,k);
        idxTest = test(c,k);
        W_test = my_ridge_regression(Y(idxTrain,:),X(idxTrain,:),lambdas(i));
        Y_test = X(idxTest,:) * W_test;
        mse_k(k) = immse(Y_test,Y(idxTest,:));
    end
    mse_i(i) = mean(mse_k);
end
[~,minIdx] = min(mse_i);
best_lambda = lambdas(minIdx);

% figure
% plot(lambdas,mse_i,'.','MarkerSize',20);
% xlabel('lambda')
% ylabel('mse')
% title([num2str(K) '-fold' ' cross-validation'])
end