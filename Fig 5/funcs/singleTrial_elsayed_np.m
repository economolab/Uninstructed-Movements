function rez = singleTrial_elsayed_np(input_data,obj,me,params,cond2use, cond2proj, nullalltime, onlyAW, delayOnly)
warning('off', 'manopt:getHessian:approx')

%% trials to use (afc)
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


%%
% trials = trialsHit; % only hit trials

if onlyAW
    trials = cell2mat(trials_cond(5:end)');
else
    trials = cell2mat(trials_cond'); % all trials from cond2use
end
%% split data into quiet and moving time points

if delayOnly
    null = [];
    potent = [];
    for trix = 1:numel(trials)
        t = trials(trix);
        delay_t(1) = obj.bp.ev.delay(t) - 2.5;
        delay_t(2) = obj.bp.ev.goCue(t) - 0.02 - 2.5;
        for i = 1:numel(delay_t)
            [~,ix(i)] = min(abs(obj.time - delay_t(i)));
        end
        
        mask{trix} = me.move(ix(1):ix(2),t);

        full{trix} = squeeze(input_data(ix(1):ix(2),t,:));
        try % if no non-move time points
            null = cat(1,null,full{trix}(~mask{trix},:));
        catch
        end
        try % if no move time points
         potent = cat(1,potent,full{trix}(mask{trix},:));
        catch
        end

    end
    N.delay = null;
    N.resp  = potent;
    N.null = null;
    N.potent = null;

    rez.N = N;

else
    % motion energy
    mask = me.move(:,trials);
    mask = mask(:); % (time*trials) , 1 where animal is moving, 0 where animal is quiet

    % single trial neural data
    N.full = input_data(:,trials,:);
    N.dims = size(N.full);
    N.full_reshape = reshape(N.full,N.dims(1)*N.dims(2),N.dims(3));

    if nullalltime
        N.null = N.full_reshape(:,:);
    else
        N.null = N.full_reshape(~mask,:);
    end

    % sample same number of null time points and potent time points
    if nullalltime
        nNull = size(N.null,1); % how many null time points
        maskix = find(mask);
        mask_ = mask(randsample(maskix,nNull,false));
    else
        mask_ = mask;
        N.potent = N.full_reshape(mask_,:);
    end


    % get delay and response epoch neural activity (only used for variance
    % explained calcs)
    delay_edges = [-0.42 -0.02];
    resp_edges  = [0.02 0.42];
    for i = 1:2
        [~,delayix(i)] = min(abs(obj(1).time - delay_edges(i)));
        [~,respix(i)] = min(abs(obj(1).time - resp_edges(i)));
    end

    N.delay = N.full(delayix(1):delayix(2),:,:);
    N.resp = N.full(respix(1):respix(2),:,:);
    N.delay = reshape(N.delay,size(N.delay,1)*size(N.delay,2),size(N.delay,3));
    N.resp = reshape(N.resp,size(N.resp,1)*size(N.resp,2),size(N.resp,3));



    rez.N = N;

end




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

rez.covDelay = cov(N.delay);
rez.covResp = cov(N.resp);

% -----------------------------------------------------------------------
% -- number of null and potent dims --
% -----------------------------------------------------------------------
% assign num dims by amount of PCs needed to explain some amount of
% variance defined in params (capped at 10 dims)
rez.varToExplain = params.N_varToExplain;
% [rez.dPrep,rez.dMove] = getNumDims(N,rez.varToExplain);
% if rez.dPrep > 10; rez.dPrep = 10; end
% if rez.dMove > 10; rez.dMove = 10; end
% if rez.dPrep == 1; rez.dPrep = 10; end
% if rez.dMove == 1; rez.dMove = 10; end
rez.dPrep = floor(size(rez.covNull,1)/2);
rez.dMove = ceil(size(rez.covNull,1)/2);
if rez.dPrep > 20; rez.dPrep = 20; end
if rez.dMove > 20; rez.dMove = 20; end

% method 2 - keep reducing var2explain until dMove+dPrep <= full dim
% check = 1;
% while check
%     [rez.dPrep,rez.dMove] = getNumDims(N,rez.varToExplain);
%     if (rez.dPrep + rez.dMove) <= size(N.full_reshape,2)
%         break
%     end
%     rez.varToExplain = rez.varToExplain - 1;
% end




% main optimization step
alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
[Q,~,P,~,~] = orthogonal_subspaces(rez.covPotent,rez.dMove, ...
    rez.covNull,rez.dPrep,alpha);


rez.Qpotent = Q*P{1};
rez.Qnull = Q*P{2};

%% projections

rez = ta_projectNP(input_data,rez,cond2proj,params);

%% var exp

rez = var_exp_NP(trials_cond,input_data,rez);

rez = var_exp_NP_recon(input_data,rez,cond2proj,params, trials, me);

end





