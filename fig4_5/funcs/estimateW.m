function rez = estimateW(obj,dat,params,time,me)

%% trials to use
% going to use hits and misses (including early), no autowater, no stim
% the later analysis will not include early trials because epochs have
% different time lengths, but here we just want to estimate these spaces

% set conditions to sue
condition(1)     = {'R&hit&~stim.enable&~autowater'};         % right hits, no stim, aw off                                                                                           enable&~autowater'};         % right hits, no stim, aw off
condition(end+1) = {'L&hit&~stim.enable&~autowater'};         % left hits, no stim, aw off
condition(end+1) = {'R&miss&~stim.enable&~autowater'};          % right hits, no stim, aw on
condition(end+1) = {'L&miss&~stim.enable&~autowater'};          % left hits, no stim, aw on


trialsArray = findTrials(obj,condition);

% use same number of each trial type

minHitTrials = cellfun(@(x) numel(x),trialsArray(1:2), 'UniformOutput',false);
minHitTrials = min(cell2mat(minHitTrials));

minMissTrials = cellfun(@(x) numel(x),trialsArray(3:4), 'UniformOutput',false);
minMissTrials = min(cell2mat(minMissTrials));


trials_hit = cellfun(@(x) randsample(x,minHitTrials), trialsArray(1:2), 'UniformOutput', false);
trialsHit = cell2mat(trials_hit);
trialsHit = trialsHit(:);

trials_miss = cellfun(@(x) randsample(x,minMissTrials), trialsArray(3:4), 'UniformOutput', false);
trialsMiss = cell2mat(trials_miss);
trialsMiss = trialsMiss(:);

trials = [trialsHit ; trialsMiss]; % (minTrials*numCond,1) vector of trials used to estimate null/potent spaces

%% time points to use 

align = mode(obj.bp.ev.(params.alignEvent));
e1 = mode(obj.bp.ev.delay) - align;
e2 = mode(obj.bp.ev.goCue) - align;
% timeix = obj.time > e1 & obj.time < e2;

timeix = logical(ones(size(obj.time)));


%% data

try
    tempN = dat.factors(:,:,trials);
catch
    'a'
end
tempV = dat.feats(:,trials,:);

% % split data into test and regression epoch 
% % test epoch is preparatory epoch
% % regression epoch is move epoch
% moveix = find(time>=params.move(1) & time<=params.move(2));
% prepix = find(time>=params.prep(1) & time<=params.prep(2));
% N_prep = N(prepix,:,:);
% N_move = N(moveix,:,:);
% V_prep = V(prepix,:,:);
% V_move = V(moveix,:,:);
% 
% % we will only use movement epoch data to estimate W
% 
% % reshape data to be (time*trials,nVars)
% N = permute(N,[1 3 2]);
% N = reshape(N,size(N,1)*size(N,2),size(N,3));
% 
% V = reshape(V,size(V,1)*size(V,2),size(V,3));
% 
% N_move = permute(N_move,[1 3 2]);
% N_move = reshape(N_move,size(N_move,1)*size(N_move,2),size(N_move,3));
% 
% V_move = reshape(V_move,size(V_move,1)*size(V_move,2),size(V_move,3));
% 
% assert(size(N_move,1)==size(V_move,1)) % N and V must have same elements in first dimension, something went wrong otherwise



% % using motion energy to label move and non move
N = permute(tempN,[1 3 2]);
N = N(timeix,:,:);
N = reshape(N,size(N,1)*numel(trials),size(N,3)); % reshape to (time*trials, neurons/factors)

V = tempV(timeix,:,:);
V = reshape(V,size(V,1)*numel(trials),size(V,3)); % reshape to (time*trials, features)

%% label time points as moving and non-moving

% move and non move times
tempme = me.data(timeix,trials);
mask = tempme(:) > (me.moveThresh);

Nnull = N(~mask,:);

Npotent = N(mask,:);

Vnull = V(~mask,:);

Vpotent = V(mask,:);


%% estimate W
% cross validate to find regularization parameter to use
% lambdas = linspace(0,10000,100);
% disp('Finding best regularization parameter, lambda, for regression')
% lambda = cross_validate(Vpotent,Npotent,lambdas);
% disp('DONE')
% lambda = 1.2112e+03; % JEB7, 4-29
% lambda = 10;
lambda = 0;

% compute transformation matrix, W, using ridge regression
W = my_ridge_regression(Vpotent,Npotent,lambda);

% for i = 1:size(Vpotent,2)
%     W(:,i) = regress(Vpotent(:,i),Npotent);
% end

tempN = dat.factors;
tempN = permute(dat.factors, [1 3 2]);
rez.N = reshape(tempN,size(tempN,1)*size(tempN,2),size(tempN,3)); % reshape to (time*trials, neurons/factors)
rez.V = reshape(dat.feats,size(dat.feats,1)*size(dat.feats,2),size(dat.feats,3)); % reshape to (time*trials, neurons/factors);
rez.W = W;
rez.lambda = lambda;
rez.moveix = mask;
rez.prepix = ~mask;



end



