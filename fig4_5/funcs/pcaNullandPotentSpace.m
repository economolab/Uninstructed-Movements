function rez = pcaNullandPotentSpace(obj,dat,me,params)

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

%% time points to use 

% align = mode(obj.bp.ev.(params.alignEvent));
% e1 = mode(obj.bp.ev.delay) - align;
% e2 = mode(obj.bp.ev.goCue) - align;
% timeix = obj.time > e1 & obj.time < e2;

timeix = logical(ones(size(obj.time)));


%% preprocess input data

% temp = dat;

temp = dat(:,:,trials); % data to estimate null/potent spaces (time,neurons/factors,trials)

N = permute(temp,[1,3,2]); % reshape to (time,trials,neurons/factors), easier to reshape with these dims


% mean center
means = zeros(size(N,1),size(N,3));
for i = 1:size(N,3)
    means(:,i) = mean(N(:,:,i),2);
end
means = repmat(means,[1 1 size(N,2)]);
means = permute(means,[1 3 2]);
N = N - means;

N = N(timeix,:,:);
N = reshape(N,size(N,1)*numel(trials),size(N,3)); % reshape to (time*trials, neurons/factors)

%% label time points as moving and non-moving

% move and non move times
tempme = me.data(timeix,trials);
mask = tempme(:) > (me.moveThresh); % when moving
% mask = reshape(mask,size(me.data,1),numel(trials));
% mask = mask(logical(timeix),:);

Nnull = N(~mask,:);

Npotent = N(mask,:);


%% null and potent spaces

rez.covNull = cov(Nnull);
rez.covPotent = cov(Npotent);

rez.varToExplain = 85;

[pcs,~,explained] = myPCA(Nnull);
rez.dPrep = numComponentsToExplainVariance(explained, rez.varToExplain );
rez.Qnull = pcs(:,1:rez.dPrep);

[~,~,explained] = myPCA(Npotent);
rez.dMove= numComponentsToExplainVariance(explained, rez.varToExplain );

% project out null space
modesToKeep = eye(size(pcs,1)) - (rez.Qnull*rez.Qnull');

proj = N * modesToKeep;

% FIND POTENT MODES 
moveproj = proj; % use all leftover data for potent mode ID (seems more right)

% pca
[pcs,~,~] = myPCA(moveproj);
rez.Qpotent = pcs(:,1:rez.dMove);



%% removing activity along the null dimensions
% https://www.sciencedirect.com/science/article/pii/S0896627319300534#sec4

% M = rez.Qnull' * cov(N);
% 
% [~,S,V] = svd(M);
% 
% % find columns of S with only 0s (corresponding cols of V make up null
% % space of M -- and M is cov of neural activity in our Null Space -- so those cols of V are our Potent Space)
% cols = find(sum(S)==0);
% 
% rez.Qpotent = V(:,cols);

%%

% project neural activity onto null and potent spaces, reshape

temp = permute(dat,[1 3 2]);
temp2 = reshape(temp,size(temp,1)*size(temp,2),size(temp,3)); % reshape to (time*trials, neurons/factors)

N_potent = temp2 * rez.Qpotent;
N_null = temp2 * rez.Qnull;

rez.N_potent = reshape(N_potent,size(temp,1),size(temp,2),min(rez.dMove,numel(cols)));
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