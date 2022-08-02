function rez = elsayedPSTH_NullandPotentSpace(obj,params,dfparams)


%% PREPROCESS
rez.psth = obj.psth(:,:,[2,3]);
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
% edges = [-0.55 -0.05];
edges = [-1.02 -0.02];
% edges = [-0.52 -0.02];
[~,e1] = min(abs(obj.time - edges(1)));
[~,e2] = min(abs(obj.time - edges(2)));
rez.prepix = e1:e2;

% move
% % edges = [.05 .55];
edges = [0.02 1.02];
% edges = [0.02 0.52];
[~,e1] = min(abs(obj.time - edges(1)));
[~,e2] = min(abs(obj.time - edges(2)));
rez.moveix = e1:e2;

%% GET COVARIANCE MATRICES of EPOCHS

% concatenates psth to (ct x n)
psthprep = rez.psth_processed(rez.prepix,:,1);
psthmove = rez.psth_processed(rez.moveix,:,1);
for i = 2:size(rez.psth_processed,3)
    psthprep = [psthprep ; rez.psth_processed(rez.prepix,:,i)];
    psthmove = [psthmove ; rez.psth_processed(rez.moveix,:,i)];
end

% computes covariance of each epoch
rez.Cprep = cov(psthprep);
rez.Cmove = cov(psthmove);

%% OPTIMIZATION

% find number of dims for move and prep epochs

[~,~,explained] = myPCA(psthprep);
rez.dPrep = numComponentsToExplainVariance(explained, 90);
if rez.dPrep == 1
    rez.dPrep = 2;
end


[~,~,explained] = myPCA(psthmove);
rez.dMove = numComponentsToExplainVariance(explained, 90);
if rez.dMove == 1
    rez.dMove = 2;
end


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

f = figure;
for i = 1:size(Q,2) % dimension
    ax = nexttile;
    title(['Potent Dim ' num2str(i)]);
%     xlim([params.tmin, params.tmax]);
    hold on
    for j = 1:2%size(obj.psth,3) % condition
        proj = rez.psth_processed(:,:,j) * Q(:,i);
        plot(obj.time, mySmooth(proj,100), 'Color', cols{j}, ...
            'LineWidth', 2.5);
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
    for j = 1:2%size(obj.psth,3) % condition
        proj = rez.psth_processed(:,:,j) * Q(:,i);
        plot(obj.time, mySmooth(proj,100), 'Color', cols{j}, ...
            'LineWidth', 2.5);
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
for i = 1:2 % condition
    rez.proj_null(:,:,i) = rez.psth(:,:,i) * Q;
end

% figure
% for i = 1:size(rez.proj_null,2)
%     ax = nexttile;
%     imagesc(squeeze(rez.proj_null(:,i,:))');
% end


Q = rez.Qpotent;

rez.proj_potent = zeros(size(rez.psth,1),size(Q,2),size(rez.psth,3)); % (time,dims,conditions)
for i = 1:2 % condition
    rez.proj_potent(:,:,i) = rez.psth(:,:,i) * Q;
end

% figure
% for i = 1:size(rez.proj_potent,2)
%     ax = nexttile;
%     imagesc(squeeze(rez.proj_potent(:,i,:))');
% end


%% single trial projections

dat{1} = obj.trialdat(:,:,params.trialid{2});
dat{2} = obj.trialdat(:,:,params.trialid{3});

projdat = cellfun(@(x) permute(x, [1 3 2]),dat,'UniformOutput',false);
projdat_reshape = cellfun(@(x) reshape(x,size(x,1)*size(x,2),size(x,3)), projdat, 'UniformOutput',false);


for i = 1:numel(dat)
    projdat_null{i} = projdat_reshape{i} * rez.Qnull;
    projdat_potent{i} = projdat_reshape{i} * rez.Qpotent;
end

for i = 1:2
    projdat_null_reshape{i} = reshape(projdat_null{i}, numel(obj.time), numel(params.trialid{i+1}), rez.dPrep);
    projdat_potent_reshape{i} = reshape(projdat_potent{i}, numel(obj.time), numel(params.trialid{i+1}), rez.dMove);
end

cols{1} = clrs.rhit;
cols{2} = clrs.lhit;
f = figure;
f.Position = [280   498   560   420];
for i = 1:rez.dPrep
    ax = nexttile; hold on;
    for j = 1:2
        for k = 1:numel(params.trialid{j+1})
            toplot = projdat_null_reshape{j}(:,k,i);
            toplot = mySmooth(toplot,31);
            patchline(obj.time,toplot,'EdgeColor',cols{j},'LineWidth',1,'EdgeAlpha',0.3)
        end
    end
    hold off
    sgtitle('Null')
end

cols{1} = clrs.rhit;
cols{2} = clrs.lhit;
f = figure;
f.Position = [844   498   560   420];
for i = 1:rez.dMove
    ax = nexttile; hold on;
    for j = 1:2
        for k = 1:numel(params.trialid{j+1})
            toplot = projdat_potent_reshape{j}(:,k,i);
            toplot = mySmooth(toplot,31);
            patchline(obj.time,toplot,'EdgeColor',cols{j},'LineWidth',1,'EdgeAlpha',0.3)       
        end
    end
    hold off
    sgtitle('Potent')
end



end

