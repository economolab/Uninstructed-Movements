function rez = get_NP_from_W(rez)
% primer on using svd to find spaces of a matrix
% http://pillowlab.princeton.edu/teaching/statneuro2018/slides/notes03a_SVDandLinSys.pdf

% rank of W (how many linearly independent cols are there)
% there will be k many potent dimensions, and size(W,2)-k null dimensions
tolerance = 0.1; % rank(A,TOL) is the number of singular values of A that are larger than TOL.
k = rank(rez.W', tolerance); 
% k = rank(rez.W);

% column, row, and null space of W can be found through SVD
[u,s,v] = svd(rez.W'); % W' = u*s*v'. check this with the command: immse(W',u*s*v'). should return ~0

% column space of W (just b/c, don't need it)
% W_col = u(:,1:k);

% row space of W
W_row = v(:,1:k);
rez.Qpotent = W_row;

% null space of W
% rez.Qnull = null(rez.W);
rez.Qnull = v(:,(k+1):end);



end 