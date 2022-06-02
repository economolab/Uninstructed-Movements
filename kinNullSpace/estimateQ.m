function rez = estimateQ(obj,dat,me,params,type)
warning('off', 'manopt:getHessian:approx')

rhit = params.trialid{2};
lhit = params.trialid{3};
rmiss = params.trialid{4};
lmiss = params.trialid{5};
trials = [rhit;lhit];

% temp = dat.(type);
temp = dat.gpfalatents;

N = permute(temp,[1,3,2]);

% mean center
means = zeros(size(N,1),size(N,3));
for i = 1:size(N,3)
    means(:,i) = mean(N(:,:,i),2);
end
means = repmat(means,[1 1 size(N,2)]);
means = permute(means,[1 3 2]);
N = N - means;
% 
% % normalize
% for i = 1:size(N,3)
%     N(:,:,i) = N(:,:,i) ./ repmat(mean(N(:,:,i),2),1,size(N,2));
% end

% N = normalize(N,2);

N_alltrials = N;
N = N(:,trials,:);

N = reshape(N,size(N,1)*numel(trials),size(N,3));


% move and non move times
tempme = me.data(:,trials);
mask = tempme(:) > (me.moveThresh);

Nnull = N(~mask,:);

Npotent = N(mask,:);


covNull = cov(Nnull);
covPotent = cov(Npotent);


% main optimization step
alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
rez.nullDims = 3;
rez.potentDims = 7;
[Q,~,P,~,~] = orthogonal_subspaces(covPotent,rez.potentDims, ...
    covNull,rez.nullDims,alpha);


rez.optim.Qpotent = Q*P{1};
rez.optim.Qnull = Q*P{2};

% project neural activity onto null and potent spaces, reshape
nTrials = obj.bp.Ntrials;
rez.N = reshape(N_alltrials,size(N_alltrials,1)*nTrials ,size(N_alltrials,3));


N_potent = rez.N * rez.optim.Qpotent;
N_null   = rez.N * rez.optim.Qnull;

N_potent = reshape(N_potent,size(temp,1),nTrials ,size(N_potent,2));
tempdat = permute(N_potent,[1,3,2]);
rez.optim.N_potent{1} = mean(tempdat(:,:,rhit),3);
rez.optim.N_potent{2} = mean(tempdat(:,:,lhit),3);

N_null = reshape(N_null,size(temp,1),nTrials ,size(N_null,2));
tempdat = permute(N_null,[1,3,2]);
rez.optim.N_null{1} = mean(tempdat(:,:,rhit),3);
rez.optim.N_null{2} = mean(tempdat(:,:,lhit),3);

rez.optim.covNull = cov(Nnull);
rez.optim.covPotent = cov(Npotent);

rez.N = N_alltrials;

% variance explained

C = cov(N); % var explained for just data used to estimate Q
nDims = size(N,2);

evals = eig(C);
evals = sort(evals,'descend');
eigsum = sum(evals(1:nDims)); 

rez.varexp_null = var_proj(rez.optim.Qnull,C,eigsum);

rez.varexp_potent = var_proj(rez.optim.Qpotent,C,eigsum);

evals = eig(rez.optim.covNull);
evals = sort(evals,'descend');
eigsum = sum(evals(1:nDims));
rez.varexp_null_nonmove = var_proj(rez.optim.Qnull,covNull,eigsum);
evals = eig(rez.optim.covPotent);
evals = sort(evals,'descend');
eigsum = sum(evals(1:nDims));
rez.varexp_null_move = var_proj(rez.optim.Qnull,covPotent,eigsum);

C = rez.optim.covPotent;
evals = eig(C);
evals = sort(evals,'descend');
eigsum = sum(evals(1:nDims));
rez.varexp_potent_move = var_proj(rez.optim.Qpotent,covPotent,eigsum);
C = rez.optim.covNull;
evals = eig(C);
evals = sort(evals,'descend');
eigsum = sum(evals(1:nDims));
rez.varexp_potent_nonmove = var_proj(rez.optim.Qpotent,covNull,eigsum);


end 