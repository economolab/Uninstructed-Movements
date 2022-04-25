function rez = calVarExp(rez)

C = cov(rez.N);
nDims = size(rez.W,1);

evals = eig(cov(rez.N));
evals = sort(evals,'descend');
eigsum = sum(evals(1:nDims));

rez.varexp_null = var_proj(rez.W_null,C,eigsum);

rez.varexp_potent = var_proj(rez.W_potent,C,eigsum);


end 