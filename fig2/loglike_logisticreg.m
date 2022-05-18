function llik = loglike_logisticreg(b, X, Y)  
num = X * b;  
prb = exp(num .* Y) ./ (1 + exp(num));  
llik = -sum(log(prb));  
end
