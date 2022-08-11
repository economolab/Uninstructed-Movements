function W = my_ridge_regression(Y,X,lambda)
% returns transformation matrix W
% inputs:
% Y - target matrix (T,M)
% X - data matrix   (T,N)
% lambda = regularization parameter

% [T,M] = size(Y);
% [~,N] = size(X);
% 
% Xplus = [X; sqrt(lambda)*eye(N)];
% Yplus = [Y; zeros(N,M)];
% 
% W = Xplus\Yplus;

W = (X'*X + lambda*eye(size(X,2))) \ (X'*Y); %(M,M) \ (N,M) = inv((M,M)) * (N,M)'


end % ridge_regression