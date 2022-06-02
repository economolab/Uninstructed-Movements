function rez = findNullPotent_pca(obj,dat,me,params,type)


rhit = params.trialid{2};
lhit = params.trialid{3};
rmiss = params.trialid{4};
lmiss = params.trialid{5};
trials = [rhit;lhit];

temp = dat.(type);

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

% null space
[pcs,~,explained] = myPCA(Nnull);
rez.dPrep = numComponentsToExplainVariance(explained, 75);
if rez.dPrep > 5
    rez.dPrep = 5;
end
rez.Qnull = pcs(:,1:rez.dPrep);

% project out null space
modesToKeep = eye(size(pcs,1)) - (rez.Qnull*rez.Qnull');

proj = N * modesToKeep;

% FIND POTENT MODES 
moveproj = proj; % use all leftover data for potent mode ID (seems more right)

% pca
[pcs,~,explained] = myPCA(moveproj);

% find number of pcs to explain 95% of variance of moveproj (not psth)
rez.dMove = numComponentsToExplainVariance(explained, 75);
if rez.dMove > 5
    rez.dMove = 5;
end
rez.Qpotent = pcs(:,1:rez.dMove);



% project neural activity onto null and potent spaces, reshape
nTrials = obj.bp.Ntrials;
rez.N = reshape(N_alltrials,size(N_alltrials,1)*nTrials ,size(N_alltrials,3));


N_potent = rez.N * rez.Qpotent;
N_null   = rez.N * rez.Qnull;

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

rez.varexp_null = var_proj(rez.Qnull,C,eigsum);

rez.varexp_potent = var_proj(rez.Qpotent,C,eigsum);

evals = eig(rez.optim.covNull);
evals = sort(evals,'descend');
eigsum = sum(evals(1:nDims));
rez.varexp_null_nonmove = var_proj(rez.Qnull,covNull,eigsum);
evals = eig(rez.optim.covPotent);
evals = sort(evals,'descend');
eigsum = sum(evals(1:nDims));
rez.varexp_null_move = var_proj(rez.Qnull,covPotent,eigsum);

C = rez.optim.covPotent;
evals = eig(C);
evals = sort(evals,'descend');
eigsum = sum(evals(1:nDims));
rez.varexp_potent_move = var_proj(rez.Qpotent,covPotent,eigsum);
C = rez.optim.covNull;
evals = eig(C);
evals = sort(evals,'descend');
eigsum = sum(evals(1:nDims));
rez.varexp_potent_nonmove = var_proj(rez.Qpotent,covNull,eigsum);



end