function ai = get_alignment_index(Q,C,sigma)

ai = trace(Q' * C * Q) ./ sum(sigma);

end