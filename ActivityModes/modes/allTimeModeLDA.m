function allmode = allTimeModeLDA(obj,cond)
% action mode: defined during mvmt init (0.1 to 0.3 s rel to go cue)
%       (hitR - hitL) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

% want a matrix X that is mean fr on each trial during delay period
% (trials,neurons)
% Y (trials,1) is trial type 

trials.num = cell(trials.N,1);
for condix = 1:trials.N
    trials.num{condix} = find(trials.ix(:,condix));
    trials.type{condix} = condix * ones(size(trials.num{condix}));
end
trialid = cell2mat(trials.num);
trialtype = cell2mat(trials.type');


X = zeros(size(obj.trialdat,2),numel(trialid));

for trix = 1:numel(trialid)
    trid = trialid(trix);
    X(:,trix) = mean(obj.trialdat(:,:,trid));
end
X = X'; % (trials,neurons)
Y = trialtype;

% calculate choice mode via lda
mLDA = LDA(X, num2cell(Y));
mLDA.Compute();

allmode = mLDA.EigenVectors(:,1);
allmode = allmode./sum(abs(allmode)); % (ncells,1)

end % actionMode