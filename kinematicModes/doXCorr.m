function mode = doXCorr(X, Y)
%% Y are observations (N x 1)
% X are predictors (N x M)
% mode is loadings (M x 1)

r = zeros(size(X, 2), 1);

for i = 1:size(X, 2)
    x = X(:, i);
    tmp = corrcoef(x, Y);
    r(i) = tmp(1,2);
end

mode = r./sum(abs(r));