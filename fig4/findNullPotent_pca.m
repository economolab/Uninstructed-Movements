function rez = findNullPotent_pca(obj,dat,me,type)


% for lfads data
rhit = find(obj.bp.R & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
lhit = find(obj.bp.L & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
mask = ismember(dat.trials,rhit);
rhit = find(mask);
mask = ismember(dat.trials,lhit);
lhit = find(mask);

trials = [rhit; lhit];
rmask = 1:numel(rhit);
lmask = (numel(rhit)+1):(numel(rhit)+1 + numel(lhit) - 1);

rez.rhit = rhit;
rez.lhit = lhit;

temp = dat.(type);

N = permute(temp,[1,3,2]);

N = N(:,trials,:);

% mean center
means = zeros(size(N,1),size(N,3));
for i = 1:size(N,3)
    means(:,i) = mean(N(:,:,i),2);
end
means = repmat(means,[1 1 size(N,2)]);
means = permute(means,[1 3 2]);
N = N - means;

tempN = reshape(N,size(N,1)*numel(trials),size(N,3));



% move and non move times
tempme = me.data(:,trials);
mask = tempme(:) > me.moveThresh;

Nnull = tempN(~mask,:);

Npotent = tempN(mask,:);


% null space
[pcs,~,explained] = myPCA(Nnull);
rez.dPrep = numComponentsToExplainVariance(explained, 95);
if rez.dPrep > 5
    rez.dPrep = 5;
end
rez.Qnull = pcs(:,1:rez.dPrep);

% project out null space
modesToKeep = eye(size(pcs,1)) - (rez.Qnull*rez.Qnull');

proj = tempN * modesToKeep;

% FIND POTENT MODES 
moveproj = proj; % use all leftover data for potent mode ID (seems more right)

% pca
[pcs,~,explained] = myPCA(moveproj);

% find number of pcs to explain 95% of variance of moveproj (not psth)
rez.dMove = numComponentsToExplainVariance(explained, 95);
if rez.dMove > 5
    rez.dMove = 5;
end
rez.Qpotent = pcs(:,1:rez.dMove);



% project neural activity onto null and potent spaces, reshape
N = reshape(N,size(N,1)*numel(trials),size(N,3));


N_potent = N * rez.Qpotent;
N_null   = N * rez.Qnull;

N_potent = reshape(N_potent,size(temp,1),numel(trials),size(N_potent,2));
rez.N_potent = permute(N_potent,[1,3,2]);

N_null = reshape(N_null,size(temp,1),numel(trials),size(N_null,2));
rez.N_null = permute(N_null,[1,3,2]);


rez.covNull = cov(Nnull);
rez.covPotent = cov(Npotent);

rez.N = N;


end