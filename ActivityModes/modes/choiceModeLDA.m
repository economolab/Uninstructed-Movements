function choicemode = choiceModeLDA(obj,params,cond,epoch,alignEvent)
% 2. choice mode: defined during delay period
%       ((hitR - missR) + (missL - hitL)) / sqrt(sum(sd for each tt ^2));

% which trials to use for each condition used for finding the mode
trials = getTrialsForModeID(obj,cond);

% find time in each trial corresponding to epoch
epochix = nan(obj.bp.Ntrials,2);
for trix = 1:obj.bp.Ntrials
    epochix(trix,:) = findedges(obj.time,obj.bp,params.dt,epoch,trix,alignEvent); % (idx1,idx2)
end

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
    X(:,trix) = mean(obj.trialdat(epochix(trix,1):epochix(trix,2),:,trid));
end
X = X'; % (trials,neurons)
Y = trialtype;

% calculate choice mode via lda
mLDA = LDA(X, num2cell(Y));
mLDA.Compute();

choicemode = mLDA.EigenVectors(:,1);
choicemode = choicemode./sum(abs(choicemode)); % (ncells,1)

end % choiceModeLDA














