function [W_null,W_potent,N_null,N_potent,dPrep,dMove] = getNullPotentSpaces_SVD(W,rez,dat)
% primer on using svd to find spaces of a matrix
% http://pillowlab.princeton.edu/teaching/statneuro2018/slides/notes03a_SVDandLinSys.pdf

% rank of W (how many linearly independent cols are there)
% there will be k many potent dimensions, and size(W,2)-k null dimensions
k = rank(W); 

% column, row, and null space of W can be found through SVD
[u,s,v] = svd(W); % W' = u*s*v'. check this with the command: immse(W',u*s*v'). should return ~0

% column space of W (just b/c, don't need it)
% W_col = u(:,1:k);

% row space of W
W_row = v(:,1:k);
W_potent = W_row;

% null space of W
W_null = null(W);

% project neural activity onto null and potent spaces, reshape
N_potent = rez.N * W_potent;
N_null   = rez.N * W_null;

N_potent = reshape(N_potent,size(dat.factors,1),size(dat.factors,3),size(N_potent,2));

N_null = reshape(N_null,size(dat.factors,1),size(dat.factors,3),size(N_null,2));

dPrep = size(N_null,3);
dMove = size(N_potent,3);



end % 