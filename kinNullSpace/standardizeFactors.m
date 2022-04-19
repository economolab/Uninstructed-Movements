function factors = standardizeFactors(factors)
factors = permute(factors,[1 3 2]);
k = reshape(factors, size(factors, 1).*size(factors, 2), size(factors, 3));
k = (k-nanmean(k, 1))./nanstd(k, [], 1);
k = reshape(k, size(factors, 1), size(factors, 2), size(factors, 3));
factors = permute(k,[1 3 2]);
end
