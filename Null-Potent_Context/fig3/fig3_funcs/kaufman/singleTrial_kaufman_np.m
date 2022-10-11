function rez = singleTrial_kaufman_np(input_data,obj,me,params,cond2use)
warning('off', 'manopt:getHessian:approx')

%% trials to use
% going to use hits and misses (including early), no autowater, no stim
% the later analysis will not include early trials because epochs have
% different time lengths, but here we just want to estimate these spaces

trials_cond = params.trialid(cond2use);

% use same number of l/r hits and same number of l/r miss

minHitTrials = cellfun(@(x) numel(x),trials_cond(1:2), 'UniformOutput',false);
nhits = min(cell2mat(minHitTrials));

minMissTrials = cellfun(@(x) numel(x),trials_cond(3:4), 'UniformOutput',false);
nmiss = min(cell2mat(minMissTrials));


trials_hit = cellfun(@(x) randsample(x,nhits), trials_cond(1:2), 'UniformOutput', false);
trialsHit = cell2mat(trials_hit);
trialsHit = trialsHit(:);

trials_miss = cellfun(@(x) randsample(x,nmiss), trials_cond(3:4), 'UniformOutput', false);
trialsMiss = cell2mat(trials_miss);
trialsMiss = trialsMiss(:);

trials = [trialsHit ; trialsMiss]; % (minTrials*numCond,1) vector of trials used to estimate null/potent spaces


% trials = trialsHit; % only hit trials


% trials = cell2mat(trials_cond'); % all trials from cond2use



%% split data into quiet and moving time points

% motion energy
mask = me.move(:,trials);
mask = mask(:); % (time*trials) , 1 where animal is moving, 0 where animal is quiet

% single trial neural data
N.full = input_data.N(:,trials,:);
N.dims = size(N.full);
N.full_reshape = reshape(N.full,N.dims(1)*N.dims(2),N.dims(3));

N.null = N.full_reshape(~mask,:);

N.potent = N.full_reshape(mask,:);

rez.N = N;

% single trial kin data
M.full = input_data.M(:,trials,:);
M.dims = size(M.full);
M.full_reshape = reshape(M.full,M.dims(1)*M.dims(2),M.dims(3));

M.null = M.full_reshape(~mask,:);

M.potent = M.full_reshape(mask,:);

rez.M = M;


%% Find W in M = WN (ridge regression)

lambda = 0.01; % should be cross-validated
rez.W = my_ridge_regression(rez.M.potent,rez.N.potent,lambda);


%% find null and potent spaces

rez = get_NP_from_W(rez);

rez.dPrep = size(rez.Qnull,2);
rez.dMove = size(rez.Qpotent,2);


%% projections

rez = projectNP(trials_cond,input_data.N,rez);

%% var exp

rez.covNull = cov(N.null);
rez.covPotent = cov(N.potent);
rez = var_exp_NP(trials_cond,input_data.N,rez);


%% trial average projections

rez = ta_projectNP(input_data.N,rez,cond2use,params);


end

% for dim = 1:10
% figure; 
% subplot(2,1,1); hold on
% plot(obj.time,squeeze(rez.N_potent_psth(:,dim,1)),'b')
% plot(obj.time,squeeze(rez.N_potent_psth(:,dim,2)),'r')
% subplot(2,1,2); hold on
% plot(obj.time,squeeze(rez.N_null_psth(:,dim,1)),'b')
% plot(obj.time,squeeze(rez.N_null_psth(:,dim,2)),'r')
% end
















