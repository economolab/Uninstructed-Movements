function [rez,times] = cdNullSpace_binnedRates(rez,obj,params,sessnum)


%% epochs

alignment = mode(obj.bp.ev.(params.alignEvent));

% early sample epoch
e1 = mode(obj.bp.ev.delay) - 0.5 - alignment; 
e2 = mode(obj.bp.ev.delay) - 0.1 - alignment;
times.early = obj.time>=e1 & obj.time <= e2;

% late delay epoch
e1 = mode(obj.bp.ev.goCue) - 0.5 - alignment; 
e2 = mode(obj.bp.ev.goCue) - 0.1 - alignment;
times.late = obj.time>=e1 & obj.time <= e2;

% go epoch    
e1 = mode(obj.bp.ev.goCue) + 0.02 - alignment; 
e2 = mode(obj.bp.ev.goCue) + 0.42 - alignment;
times.go = obj.time>=e1 & obj.time <= e2;

%% psths (trial-averaged null and potent space projections)

rez.null_psth = zeros(size(rez.N_null,1),size(rez.N_null,3),2); % (time,dPrep,numCond)
rez.null_psth(:,:,1) = squeeze(mean(rez.N_null(:,params.trialid{sessnum}{2},:),2));
rez.null_psth(:,:,2) = squeeze(mean(rez.N_null(:,params.trialid{sessnum}{3},:),2));

rez.potent_psth = zeros(size(rez.N_potent,1),size(rez.N_potent,3),2); % (time,dPrep,numCond)
rez.potent_psth(:,:,1) = squeeze(mean(rez.N_potent(:,params.trialid{sessnum}{2},:),2));
rez.potent_psth(:,:,2) = squeeze(mean(rez.N_potent(:,params.trialid{sessnum}{3},:),2));

%% cd early

% null
tempdat = rez.null_psth;
mu = squeeze(mean(tempdat(times.early,:,:),1));
sd = squeeze(std(tempdat(times.early,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.null.cdEarly_mode = cd;

% potent
tempdat = rez.potent_psth;
mu = squeeze(mean(tempdat(times.early,:,:),1));
sd = squeeze(std(tempdat(times.early,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.potent.cdEarly_mode = cd;

%% cd late

% null
tempdat = rez.null_psth;
mu = squeeze(mean(tempdat(times.late,:,:),1));
sd = squeeze(std(tempdat(times.late,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.null.cdLate_mode = cd;

% potent
tempdat = rez.potent_psth;
mu = squeeze(mean(tempdat(times.late,:,:),1));
sd = squeeze(std(tempdat(times.late,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.potent.cdLate_mode = cd;

%% cd go

% null
tempdat = rez.null_psth;
mu = squeeze(mean(tempdat(times.go,:,:),1));
sd = squeeze(std(tempdat(times.go,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.null.cdGo_mode = cd;

% potent
tempdat = rez.potent_psth;
mu = squeeze(mean(tempdat(times.go,:,:),1));
sd = squeeze(std(tempdat(times.go,:,:),[],1));
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.potent.cdGo_mode = cd;

%% orthogonalize

try
    rez.Qpotent = rez.W_potent;
    rez.Qnull = rez.W_null;
catch
end

[fns,~] = patternMatchCellArray(fieldnames(rez.cd.null),{'mode'},'all');
nullmodes = zeros(rez.dPrep,numel(fns));
potentmodes = zeros(size(rez.Qpotent,2),numel(fns));
for i = 1:numel(fns)
    nullmodes(:,i) = rez.cd.null.(fns{i});
    potentmodes(:,i) = rez.cd.potent.(fns{i});
end

orthNullModes = gschmidt(nullmodes);
orthPotentModes = gschmidt(potentmodes);

for i = 1:numel(fns)
    rez.cd.null.(fns{i}) = orthNullModes(:,i);
    rez.cd.potent.(fns{i}) = orthPotentModes(:,i);
end


%% projections

temp_potent = rez.N_potent;
temp_potent_reshape = reshape(temp_potent,size(temp_potent,1)*size(temp_potent,2),size(temp_potent,3));

temp_null = rez.N_null;
temp_null_reshape = reshape(temp_null,size(temp_null,1)*size(temp_null,2),size(temp_null,3));

for i = 1:numel(fns)
    proj_null = temp_null_reshape * rez.cd.null.(fns{i});
    rez.cd.null.([fns{i}(1:end-5) '_latent']) = reshape(proj_null,size(temp_null,1),size(temp_null,2),1);
    
    proj_potent = temp_potent_reshape * rez.cd.potent.(fns{i});
    rez.cd.potent.([fns{i}(1:end-5) '_latent']) = reshape(proj_potent,size(temp_potent,1),size(temp_potent,2),1);
end


%% var exp

temp = cat(1,rez.null_psth(:,:,1),rez.null_psth(:,:,2));
nullcov = cov(temp);
[~,eigvals] = myPCA(nullcov);
% eigvals = 1;
rez.cd.null.ve.early = var_proj(rez.cd.null.cdEarly_mode,nullcov,sum(eigvals));
rez.cd.null.ve.late = var_proj(rez.cd.null.cdLate_mode,nullcov,sum(eigvals));
rez.cd.null.ve.go = var_proj(rez.cd.null.cdGo_mode,nullcov,sum(eigvals));

temp = cat(1,rez.potent_psth(:,:,1),rez.potent_psth(:,:,2));
potentcov = cov(temp);
[~,eigvals] = myPCA(potentcov);
% eigvals = 1;
rez.cd.potent.ve.early = var_proj(rez.cd.potent.cdEarly_mode,potentcov,sum(eigvals));
rez.cd.potent.ve.late = var_proj(rez.cd.potent.cdLate_mode,potentcov,sum(eigvals));
rez.cd.potent.ve.go = var_proj(rez.cd.potent.cdGo_mode,potentcov,sum(eigvals));

end

















