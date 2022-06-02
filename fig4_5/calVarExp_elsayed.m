function rez = calVarExp_elsayed(rez)

C = cov(rez.N);
nDims = size(rez.N,2);

evals = eig(cov(rez.N));
evals = sort(evals,'descend');
eigsum = sum(evals(1:nDims));

rez.varexp_null = var_proj(rez.optim.Qnull,C,eigsum);

rez.varexp_potent = var_proj(rez.optim.Qpotent,C,eigsum);


end 