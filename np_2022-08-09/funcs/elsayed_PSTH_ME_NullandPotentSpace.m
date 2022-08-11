function rez = elsayed_PSTH_ME_NullandPotentSpace(obj,me,params,dfparams,kin,kinfeats)

temp = me.data';
temp(isnan(temp)) = 0;

dat{1} = temp(params.trialid{2},:); % right trials (trials,time)
dat{2} = temp(params.trialid{3},:); % left trials  (trials,time)

% times = [-1.0 -0.01];
% times = [-1.7 -0.01];
% times = [-2.5 -0.01];
times = [-1 -0.001];
[~,ix1] = min(abs(obj.time - times(1)));
[~,ix2] = min(abs(obj.time - times(2)));

cludat = cellfun(@(x) x(:,ix1:ix2),dat,'UniformOutput',false);

k = 5; % number of clusters;

idx = cellfun(@(x) kmeans(x,k,'Replicates',1000), cludat, 'UniformOutput',false);

tt = 1;
[sorted{tt},sortix{tt}] = sort(idx{tt});
dat_sorted{tt} = dat{tt}(sortix{tt},:);

tt = 2;
[sorted{tt},sortix{tt}] = sort(idx{tt});
dat_sorted{tt} = dat{tt}(sortix{tt},:);

f = figure;
f.Position = [263         191        1183         777];
for i = 1:numel(dat_sorted)
    subplot(1,2,i); hold on;
    tempdat = dat_sorted{i};
    tempdat = tempdat - me.moveThresh;
    tempdat(tempdat<=0) = -max(max(tempdat));
    imagesc(obj.time,1:size(dat_sorted{i}), tempdat)
%     imagesc(obj.time,1:size(dat_sorted{i}), dat_sorted{i})
    colorbar
    for j = 1:(k-1)
        trix = find(sorted{i} == j,1,'last');
        yline(trix,'w--','LineWidth',2)
    end
end

rez = 1;
return



figure;
for i = 1:2
    subplot(1,2,i)
    histogram(dat_sorted{i})
end


% me_avg = zeros(size(me.data,1),k); % (time,clusters)
% for i = 1:k
%     trix = find(idx == i);
%     me_avg(:,i) = mySmooth(nanmean(me.data(:,trix),2),31);
% end
% 
% figure; plot(obj.time,me_avg,'LineWidth',2)




%% older version of previous section
% temp = me.data';
% temp(isnan(temp)) = 0;
% 
% jawpos = kinfeats(:,:,14)'; % jaw y pos view 1
% jawpos = fillmissing(jawpos,'nearest');
% 
% % times = [-1.0 -0.01];
% times = [-1.7 -0.01];
% % times = [-2.5 -0.01];
% % times = [-2.5 2.5];
% [~,ix1] = min(abs(obj.time - times(1)));
% [~,ix2] = min(abs(obj.time - times(2)));
% 
% clutemp = temp(:,ix1:ix2);
% % clutemp = cat(2,clutemp,jawpos(:,ix1:ix2));
% 
% 
% k = 6; % number of clusters;
% 
% idx = kmeans(clutemp,k,'Replicates',100);
% 
% [sorted,sortix] = sort(idx);
% 
% temp_sorted = temp(sortix,:);
% 
% figure; imagesc(obj.time,1:obj.bp.Ntrials,temp_sorted);
% hold on;
% for i = 1:(k-1)
%     trix = find(sorted == i,1,'last');
%     yline(trix,'w--','LineWidth',2);
% end
% 
% 
% me_avg = zeros(size(me.data,1),k); % (time,clusters)
% for i = 1:k
%     trix = find(idx == i);
%     me_avg(:,i) = mySmooth(nanmean(me.data(:,trix),2),31);
% end
% 
% figure; plot(obj.time,me_avg,'LineWidth',2)


%% PREPROCESS

for i = 1:k
    trix = find(idx == i);
    temp = obj.trialdat(:,:,trix);

    rez.psth(:,:,i) = mySmooth(mean(temp,3),21);
end

rez.time = obj.time;


% soft-normalize 
rez.softnorm_lambda = 0.01;
% firing rate of a neuron, x of size (time,1), is transformed as:
% x_norm = x / (lambda + max(x) - min(x))
snpsth = rez.psth ./ (rez.softnorm_lambda + max(rez.psth) - min(rez.psth));

% mean center across conditions
for clu = 1:size(snpsth,2)
    % find mean at each time point for each condition (time,1)
    mean_across_cond = mean(snpsth(:,clu,:),3);
    snpsth(:,clu,:) = snpsth(:,clu,:) - repmat(mean_across_cond,[1,1,size(snpsth,3)]);
end
rez.psth_processed = snpsth; % soft normed and mean centered

%% PREP and MOVE EPOCHS
% prepix and moveix corresponds to time idx for each epoch

% prep
% edges = [-0.55 .05];
edges = [-2.2 0.01];
[~,e1] = min(abs(obj.time - edges(1)));
[~,e2] = min(abs(obj.time - edges(2)));
rez.prepix = e1:e2;

% move
% edges = [.05 .55];
edges = [0.01 2.2];
[~,e1] = min(abs(obj.time - edges(1)));
[~,e2] = min(abs(obj.time - edges(2)));
rez.moveix = e1:e2;

%% GET COVARIANCE MATRICES of EPOCHS

% concatenates psth to (ct x n)
psthprep = rez.psth_processed(rez.prepix,:,1);
psthmove = rez.psth_processed(rez.moveix,:,1);
for i = 2:size(rez.psth,3)
    psthprep = [psthprep ; rez.psth_processed(rez.prepix,:,i)];
    psthmove = [psthmove ; rez.psth_processed(rez.moveix,:,i)];
end

% computes covariance of each epoch
rez.Cprep = cov(psthprep);
rez.Cmove = cov(psthmove);

%% OPTIMIZATION

% find number of dims for move and prep epochs

[~,~,explained] = myPCA(psthprep);
rez.dPrep = numComponentsToExplainVariance(explained, 75);

[~,~,explained] = myPCA(psthmove);
rez.dMove = numComponentsToExplainVariance(explained, 75);

rez.dMax = max(rez.dMove,rez.dPrep);

% main optimization step
rez.alpha = 0; % regularization hyperparam (+ve->discourage sparity, -ve->encourage sparsity)
[Q,~,P,~,~] = orthogonal_subspaces(rez.Cmove,rez.dMove, ... 
                                    rez.Cprep,rez.dPrep,rez.alpha);
                                

rez.Qpotent = Q*P{1};
rez.Qnull = Q*P{2};


%% plots


Q = rez.Qpotent;

clrs = getColors;
cols{1} = clrs.rhit;
cols{2} = clrs.lhit;
cols{3} = [1, 140, 52] ./ 255;

f = figure;
for i = 1:size(Q,2) % dimension
    ax = nexttile;
    title(['Potent Dim ' num2str(i)]);
%     xlim([params.tmin, params.tmax]);
    hold on
    for j = 1:k%size(obj.psth,3) % condition
        proj = rez.psth_processed(:,:,j) * Q(:,i);
        plot(obj.time, mySmooth(proj,100), 'LineWidth', 2.5);
    end
    hold off
    ax.FontSize = 12;
end




Q = rez.Qnull;

f = figure;
for i = 1:size(Q,2) % dimension
    ax = nexttile;
    title(['Null Dim ' num2str(i)]);
%     xlim([params.tmin, params.tmax]);
    hold on
    for j = 1:k%size(obj.psth,3) % condition
        proj = rez.psth_processed(:,:,j) * Q(:,i);
        plot(obj.time, mySmooth(proj,100),'LineWidth', 2.5);
    end
    hold off
    ax.FontSize = 12;
end


%% ve
prepeigs = sort(eig(rez.Cprep),'descend');
moveeigs = sort(eig(rez.Cmove),'descend');

prepproj = rez.Qnull'*rez.Cprep*rez.Qnull;
moveproj = rez.Qpotent'*rez.Cmove*rez.Qpotent;

crossprepproj = rez.Qpotent'*rez.Cprep*rez.Qpotent;
crossmoveproj = rez.Qnull'*rez.Cmove*rez.Qnull;

rez.Qprep_null_ve = trace(prepproj) / sum(prepeigs(1:rez.dPrep)) * 100;
rez.Qmove_potent_ve = trace(moveproj) / sum(moveeigs(1:rez.dMove)) * 100;

rez.Qprep_potent_ve = trace(crossprepproj) / sum(prepeigs(1:rez.dPrep)) * 100;
rez.Qmove_null_ve = trace(crossmoveproj) / sum(moveeigs(1:rez.dMove)) * 100;

figure; 
X = categorical({'PrepNull','PrepPotent','MovePotent','MoveNull'});
Y = [rez.Qprep_null_ve,rez.Qprep_potent_ve,rez.Qmove_potent_ve,rez.Qmove_null_ve]';
bar(X,Y)

%% projections

Q = rez.Qnull;

rez.proj_null = zeros(size(rez.psth,1),size(Q,2),size(rez.psth,3)); % (time,dims,conditions)
for i = 1:k % condition
    rez.proj_null(:,:,i) = rez.psth(:,:,i) * Q;
end

% figure
% for i = 1:size(rez.proj_null,2)
%     ax = nexttile;
%     imagesc(squeeze(rez.proj_null(:,i,:))');
% end


Q = rez.Qpotent;

rez.proj_potent = zeros(size(rez.psth,1),size(Q,2),size(rez.psth,3)); % (time,dims,conditions)
for i = 1:k % condition
    rez.proj_potent(:,:,i) = rez.psth(:,:,i) * Q;
end

% figure
% for i = 1:size(rez.proj_potent,2)
%     ax = nexttile;
%     imagesc(squeeze(rez.proj_potent(:,i,:))');
% end





end

