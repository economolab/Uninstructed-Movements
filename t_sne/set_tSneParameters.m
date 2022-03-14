function parameters = set_tSneParameters(parameters)
%setRunParameters sets all parameters for the algorithms used here.
%       Any parameters not explicitly set will revert to their listed
%       default values.
%
%
% (C) Gordon J. Berman, 2014
%     Princeton University


if nargin < 1
    parameters = [];
end



%%%%%%%% t-SNE Parameters %%%%%%%%


%2^H (H is the transition entropy)
parameters.perplexity = 32;

%relative convergence criterium for t-SNE
parameters.relTol = 1e-4;

%number of dimensions for use in t-SNE
parameters.num_tsne_dim = 2;

%binary search tolerance for finding pointwise transition region
parameters.sigmaTolerance = 1e-5;

%maximum number of non-zero neighbors in P
parameters.maxNeighbors = 200;

%initial momentum
parameters.momentum = .5;

%value to which momentum is changed
parameters.final_momentum = 0.8;

%iteration at which momentum is changed
parameters.mom_switch_iter = 250;

%iteration at which lying about P-values is stopped
parameters.stop_lying_iter = 125;

%degree of P-value expansion at early iterations
parameters.lie_multiplier = 4;

%maximum number of iterations
parameters.max_iter = 1000;

%initial learning rate
parameters.epsilon = 500;

%minimum gain for delta-bar-delta
parameters.min_gain = .01;

%readout variable for t-SNE
parameters.tsne_readout = 1;

%embedding batchsize
parameters.embedding_batchSize = 20000;

%maximum number of iterations for the Nelder-Mead algorithm
parameters.maxOptimIter = 100;

%number of points in the training set
parameters.trainingSetSize = 35000;

%local neighborhood definition in training set creation
parameters.kdNeighbors = 5;

%t-SNE training set stopping critereon
parameters.training_relTol = 2e-3;

%t-SNE training set perplexity
parameters.training_perplexity = 20;

%number of points to evaluate in each training set file
parameters.training_numPoints = 10000;

%minimum training set template length
parameters.minTemplateLength = 1;
























