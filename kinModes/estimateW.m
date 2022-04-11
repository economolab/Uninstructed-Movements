function [W,N,V] = estimateW(dat,params)

switch params.full_or_reduced
    case 'full'
        N = dat.rates;
        V = dat.feats;
    case 'reduced'
        N = dat.factors;
        V = dat.feats_reduced;
    otherwise
        error('params.full_or_reduced must be set to either `full` or `reduced`')
end

% reshape data to be (time*trials,nVars)
N = permute(N,[1 3 2]);
N = reshape(N,size(N,1)*size(N,2),size(N,3));

V = permute(V,[1 3 2]);
V = reshape(V,size(V,1)*size(V,2),size(V,3));

assert(size(N,1)==size(V,1)) % N and V must have same elements in first dimension, something went wrong otherwise


% standardize neural and video data
N = zscore(N);
V = zscore(V);


% cross validate to find regularization parameter to use
lambdas = linspace(0,10000,1000);
disp('Finding best regularization parameter, lambda, for regression')
lambda = cross_validate(N,V,lambdas);
disp('DONE')
% lambda = 1.2112e+03; % JEB7, 4-29
% lambda = 10;

% compute transformation matrix, W, using ridge regression
W = my_ridge_regression(V,N,lambda);




end



