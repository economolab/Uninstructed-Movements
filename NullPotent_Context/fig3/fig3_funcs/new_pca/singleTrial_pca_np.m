function rez = singleTrial_pca_np(input_data,obj,me,params,cond2use,first)
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
N.full = input_data(:,trials,:);
N.dims = size(N.full);
N.full_reshape = reshape(N.full,N.dims(1)*N.dims(2),N.dims(3));

N.null = N.full_reshape(~mask,:);

N.potent = N.full_reshape(mask,:);


rez.N = N;


%% null and potent spaces

% -----------------------------------------------------------------------
% -- compute covariance matrices --
% -----------------------------------------------------------------------

% % method 1 - recover covariance estimate from factor analysis
% [lambda,psi] = factoran(N.null,10);
% rez.covNull = lambda*lambda' + psi;
% [lambda,psi] = factoran(N.potent,10);
% rez.covPotent = lambda*lambda' + psi;

% % method 2 - standard method
rez.covNull = cov(N.null);
rez.covPotent = cov(N.potent);

% -----------------------------------------------------------------------
% -- number of null and potent dims --
% -----------------------------------------------------------------------
% assign num dims by amount of PCs needed to explain some amount of
% variance defined in params (capped at 10 dims)
rez.varToExplain = params.N_varToExplain;
[rez.dPrep,rez.dMove] = getNumDims(N,rez.varToExplain);
if rez.dPrep > 10; rez.dPrep = 10; end
if rez.dMove > 10; rez.dMove = 10; end
if rez.dPrep == 1; rez.dPrep = 10; end
if rez.dMove == 1; rez.dMove = 10; end

% method 2 - keep reducing var2explain until dMove+dPrep <= full dim
% check = 1;
% while check
%     [rez.dPrep,rez.dMove] = getNumDims(N,rez.varToExplain);
%     if (rez.dPrep + rez.dMove) <= size(N.full_reshape,2)
%         break
%     end
%     rez.varToExplain = rez.varToExplain - 1;
% end



%% find null and potent spaces

if strcmpi(first,'null')
    second = 'potent';
    ndims1 = rez.dPrep;
    ndims2 = rez.dMove;
else
    second = 'null';
    ndims1 = rez.dMove;
    ndims2 = rez.dPrep;
end

% find first space
rez.(['Q' first]) = pca(rez.N.(first),'NumComponents',ndims1);

% find second space

% project out first space from full data
modesToKeep = eye(size(rez.(['Q' first]),1)) - (rez.(['Q' first])*rez.(['Q' first])');
proj_out = rez.N.full_reshape * modesToKeep;

% potent space
rez.(['Q' second]) = pca(proj_out,'NumComponents',ndims2);


%% projections

rez = projectNP(trials_cond,input_data,rez);

%% var exp

rez = var_exp_NP(trials_cond,input_data,rez);


%% trial average projections

rez = ta_projectNP(input_data,rez,cond2use,params);


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
















