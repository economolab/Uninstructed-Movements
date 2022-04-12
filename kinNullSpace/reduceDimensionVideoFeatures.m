function feats_reduced = reduceDimensionVideoFeatures(feats,varToExplain,nNeuralDims)

feats_rsz = permute(feats, [1 3 2]);                          % Re-order the dimensions (switch the 2nd and 3rd dimensions)
feats_rsz = reshape(feats_rsz, size(feats_rsz, 1)*size(feats_rsz, 2), size(feats_rsz, 3));   % (time*trials x feats)

% how many PCs needed to explain 90% variance
feats_rsz(isnan(feats_rsz)) = 0;
[~, ~, ~, ~, explained] = pca(feats_rsz);
numFactors = numComponentsToExplainVariance(explained, varToExplain);
% but limit numFactors to nNeuralDims - 3 (this ensures at least 3 null
% dimensions will be identified
if numFactors>(nNeuralDims-3)
    numFactors = nNeuralDims-3;
end

% [~,feats_reduced]  = pca(feats_rsz,'NumComponents',numFactors);       % Reduce the number of dimensions

[~, ~, ~, ~, feats_reduced] = factoran(feats_rsz, numFactors);
feats_reduced = reshape(feats_reduced,size(feats,1),size(feats,3),size(feats_reduced,2));
feats_reduced = permute(feats_reduced,[1,3,2]); % (time,factors,trials)

end