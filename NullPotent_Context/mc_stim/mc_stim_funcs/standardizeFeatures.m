function feats = standardizeFeatures(feats)
k = reshape(feats, size(feats, 1).*size(feats, 2), size(feats, 3));
k = (k-nanmean(k, 1))./nanstd(k, [], 1);
k = reshape(k, size(feats, 1), size(feats, 2), size(feats, 3));
feats = k;
end
