function rez = elsayedNullandPotentSpace(obj,dat,me,params)
warning('off', 'manopt:getHessian:approx')

%% trials to use
% going to use hits and misses (including early), no autowater, no stim
% the later analysis will not include early trials because epochs have
% different time lengths, but here we just want to estimate these spaces

% set conditions to sue
condition(1)     = {'R&hit&~stim.enable&~autowater'};         % right hits, no stim, aw off
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
% trials = trialsHit;

%% time points to use 

timeix = logical(ones(size(obj.time)));


%% preprocess input data

% temp = dat;

temp = dat(:,:,trials); % data to estimate null/potent spaces (time,neurons/factors,trials)

N = permute(temp,[1,3,2]); % reshape to (time,trials,neurons/factors), easier to reshape with these dims

% % for each neuron, subtract mean across all trials
% for cluix = 1:size(N,3)
%     cludat = N(:,:,cluix);
%     means = mean(cludat,2);
%     cludat = cludat - means;
%     N(:,:,cluix) = cludat;
% end


N = N(timeix,:,:);
N = reshape(N,size(N,1)*numel(trials),size(N,3)); % reshape to (time*trials, neurons/factors)

N = zscore(mySmooth(N,31));

% figure; imagesc(N')
% colorbar
% caxis([0 5])


%% label time points as moving and non-moving

% use saved movement indices
mask = me.moveIx(timeix,trials);
mask = mask(:);

% % single trials prep and move epochs
% [~,e1] = min(abs(obj.time - 0));
% [~,e2] = min(abs(obj.time - (-1.5)));
% [~,e3] = min(abs(obj.time - 1));
% me.moveIx(1:e1-1,:) = 0;
% me.moveIx(e1:end,:) = 1;
% mask = me.moveIx(timeix,trials);
% mask = mask(:);

% % N = (time,trials,factors/neurons)
% % only use +-50 ms surrounding go cue
% [~,e1] = min(abs(obj.time - (-0.05)));
% [~,e2] = min(abs(obj.time - 0.05));
% 
% N = N(e1:e2,:,:);
% N = reshape(N,size(N,1)*numel(trials),size(N,3)); % reshape to (time*trials, neurons/factors)
% mask = mask(e1:e2,:);
% mask = mask(:);


% % use move threshold and motion energy to label movment indices
% tempme = me.data(timeix,trials);
% mask = tempme(:) > (me.moveThresh);

% % shift me relative to neural data (we thougt there might have been an
% % offset but there doesn't seem to be)
% k = 6;
% newmask = [false(k, 1); mask];
% newmask = newmask(1:end-k);
% Nnull = N(~newmask,:);
% Npotent = N(newmask,:);



Nnull = N(~mask,:);

Npotent = N(mask,:);

%% null and potent spaces

rez.covNull = cov(Nnull);
rez.covPotent = cov(Npotent);

rez.varToExplain = 75;

[~,~,explained] = myPCA(Nnull);
rez.dPrep = numComponentsToExplainVariance(explained, rez.varToExplain );

[~,~,explained] = myPCA(Npotent);
rez.dMove= numComponentsToExplainVariance(explained, rez.varToExplain );

% main optimization step
alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
[Q,~,P,~,~] = orthogonal_subspaces(rez.covPotent,rez.dMove, ...
                                   rez.covNull,rez.dPrep,alpha);


rez.Qpotent = Q*P{1};
rez.Qnull = Q*P{2};

%%

% project neural activity onto null and potent spaces, reshape

temp = permute(dat,[1 3 2]);
temp2 = reshape(temp,size(temp,1)*size(temp,2),size(temp,3)); % reshape to (time*trials, neurons/factors)

N_potent = temp2 * rez.Qpotent;
N_null = temp2 * rez.Qnull;

rez.N_potent = reshape(N_potent,size(temp,1),size(temp,2),rez.dMove);
rez.N_null = reshape(N_null,size(temp,1),size(temp,2),rez.dPrep);

%% 
% project mean-centered data and store these separately. Just using this to
% see if we get qualitively different projs when projecting mean-centered
% data vs. not mean-centered (it looks like there's not much difference)


% N_reshape = reshape(N,size(temp,1),numel(trials),size(N,2));
% N_rhit = N_reshape(:,1:numel(trials_hit{1}),:);
% N_rhit = reshape(N_rhit,size(N_rhit,1)*size(N_rhit,2),size(N_rhit,3));
% rez.N_rhit_potent = N_rhit * rez.Qpotent;
% rez.N_rhit_null = N_rhit * rez.Qnull;
% rez.N_rhit_potent = reshape(rez.N_rhit_potent,size(temp,1),minHitTrials,rez.dMove);
% rez.N_rhit_null = reshape(rez.N_rhit_null,size(temp,1),minHitTrials,rez.dPrep);
% 
% N_lhit = N_reshape(:,(numel(trials_hit{1})+1):(numel(cell2mat(trials_hit))),:);
% N_lhit = reshape(N_lhit,size(N_lhit,1)*size(N_lhit,2),size(N_lhit,3));
% rez.N_lhit_potent = N_lhit * rez.Qpotent;
% rez.N_lhit_null = N_lhit * rez.Qnull;
% rez.N_lhit_potent = reshape(rez.N_lhit_potent,size(temp,1),minHitTrials,rez.dMove);
% rez.N_lhit_null = reshape(rez.N_lhit_null,size(temp,1),minHitTrials,rez.dPrep);

%% var exp

% var exp of all data used to estimate spaces
[~,eigvals,~] = myPCA(N);
rez.ve.null_total = var_proj(rez.Qnull,cov(N),sum(eigvals));
rez.ve.potent_total = var_proj(rez.Qpotent,cov(N),sum(eigvals));

% var exp of prep activity in null space
[~,eigvals,~] = myPCA(Nnull);
rez.ve.null_prep = var_proj(rez.Qnull,rez.covNull,sum(eigvals(1:rez.dPrep)));

% var exp of move activity in null space
[~,eigvals,~] = myPCA(Npotent);
rez.ve.null_move = var_proj(rez.Qnull,rez.covPotent,sum(eigvals(1:rez.dPrep)));

% var exp of move activity in potent space
[~,eigvals,~] = myPCA(Npotent);
rez.ve.potent_move = var_proj(rez.Qpotent,rez.covPotent,sum(eigvals(1:rez.dMove)));

% var exp of prep activity in potent space
[~,eigvals,~] = myPCA(Npotent);
rez.ve.potent_prep = var_proj(rez.Qpotent,rez.covNull,sum(eigvals(1:rez.dMove)));



end 