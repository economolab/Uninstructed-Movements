function rez = estimateW(dat,params,time)

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

% split data into test and regression epoch 
% test epoch is preparatory epoch
% regression epoch is move epoch
moveix = find(time>=params.move(1) & time<=params.move(2));
prepix = find(time>=params.prep(1) & time<=params.prep(2));
N_prep = N(prepix,:,:);
N_move = N(moveix,:,:);
V_prep = V(prepix,:,:);
V_move = V(moveix,:,:);

% we will only use movement epoch data to estimate W

% reshape data to be (time*trials,nVars)
N = permute(N,[1 3 2]);
N = reshape(N,size(N,1)*size(N,2),size(N,3));

V = reshape(V,size(V,1)*size(V,2),size(V,3));

N_move = permute(N_move,[1 3 2]);
N_move = reshape(N_move,size(N_move,1)*size(N_move,2),size(N_move,3));

V_move = reshape(V_move,size(V_move,1)*size(V_move,2),size(V_move,3));

assert(size(N_move,1)==size(V_move,1)) % N and V must have same elements in first dimension, something went wrong otherwise


% cross validate to find regularization parameter to use
lambdas = linspace(0,10000,1000);
disp('Finding best regularization parameter, lambda, for regression')
% lambda = cross_validate(V_move,N_move,lambdas);
% disp('DONE')
% lambda = 1.2112e+03; % JEB7, 4-29
% lambda = 10;
lambda = 0;

% compute transformation matrix, W, using ridge regression
W = my_ridge_regression(V_move,N_move,lambda);

rez.N = N;
rez.V = V;
rez.W = W;
rez.lambda = lambda;
rez.moveix = moveix;
rez.prepix = prepix;



end



