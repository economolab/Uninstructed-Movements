function rez = estimateQ(obj,dat,rez,me)


% for lfads data
rhit = find(obj.bp.R & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
lhit = find(obj.bp.L & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
mask = ismember(dat.trials,rhit);
rhit = find(mask);
mask = ismember(dat.trials,lhit);
lhit = find(mask);

trials = [rhit; lhit];
rmask = 1:numel(rhit);
lmask = (numel(rhit)+1):(numel(rhit)+1 + numel(lhit) - 1);

N = permute(dat.factors,[1,3,2]);

N = N(:,trials,:);

% mean center
means = zeros(size(N,1),size(N,3));
for i = 1:size(N,3)
    means(:,i) = mean(N(:,:,i),2);
end
means = repmat(means,[1 1 size(N,2)]);
means = permute(means,[1 3 2]);
N = N - means;

tempN = reshape(N,size(N,1)*numel(trials),size(N,3));


% move and non move times
tempme = me.data(:,trials);
mask = tempme(:) > me.moveThresh;

Nnull = tempN(~mask,:);

Npotent = tempN(mask,:);


covNull = cov(Nnull);
covPotent = cov(Npotent);


% main optimization step
alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
nullDims = 3;
potentDims = 3;
[Q,~,P,~,~] = orthogonal_subspaces(covPotent,potentDims, ...
    covNull,nullDims,alpha);


rez.optim.Qpotent = Q*P{1};
rez.optim.Qnull = Q*P{2};

% project neural activity onto null and potent spaces, reshape
N = reshape(N,size(N,1)*numel(trials),size(N,3));


N_potent = N * rez.optim.Qpotent;
N_null   = N * rez.optim.Qnull;

N_potent = reshape(N_potent,size(dat.factors,1),numel(trials),size(N_potent,2));
rez.optim.N_potent = permute(N_potent,[1,3,2]);

N_null = reshape(N_null,size(dat.factors,1),numel(trials),size(N_null,2));
rez.optim.N_null = permute(N_null,[1,3,2]);


rez.optim.covNull = cov(Nnull);
rez.optim.covPotent = cov(Npotent);




end 