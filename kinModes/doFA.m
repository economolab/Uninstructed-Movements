function dat = doFA(obj,params)

trials = cellfun(@(x)x(:),params.trialid,'un',0);
dat.trials = cat(1,trials{:});

fr = obj.trialdat(:,:,dat.trials);

fr = permute(fr, [ 1 3 2]);                          % Re-order the dimensions of filtfr (switch the 2nd and 3rd dimensions)
fr = reshape(fr, size(fr, 1)*size(fr, 2), size(fr, 3));   % (time*trials x cells)

% how many factors needed to explain 90% variance
[~, ~, ~, ~, explained] = pca(fr);
numFactors = numComponentsToExplainVariance(explained, 75);
% but limit numFactors to 10
if numFactors>10
    numFactors = 10;
end

[~,~,~,~,dat.factors]  = factoran(fr,numFactors);       % Reduce the number of neural dimensions

dat.rates = obj.trialdat(:,:,dat.trials);

end