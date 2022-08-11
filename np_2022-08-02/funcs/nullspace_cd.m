function [rez,times] = nullspace_cd(rez,obj,params)

switching = 0;
orthog = 1;

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
rez.null_psth(:,:,1) = squeeze(mean(rez.N_null(:,params.trialid{2},:),2));
rez.null_psth(:,:,2) = squeeze(mean(rez.N_null(:,params.trialid{3},:),2));

rez.potent_psth = zeros(size(rez.N_potent,1),size(rez.N_potent,3),2); % (time,dPrep,numCond)
rez.potent_psth(:,:,1) = squeeze(mean(rez.N_potent(:,params.trialid{2},:),2));
rez.potent_psth(:,:,2) = squeeze(mean(rez.N_potent(:,params.trialid{3},:),2));

%% cd early

% null
tempdat = rez.null_psth;
mu = squeeze(mean(tempdat(times.early,:,:),1));

if switching
    toswitch = find(mu(:,1) < mu(:,2));
    for i = 1:numel(toswitch)
        new = flip(mu(toswitch(i),:));
        mu(toswitch(i),:) = new;
    end
end

sd = squeeze(std(tempdat(times.early,:,:),[],1));
if switching
    for i = 1:numel(toswitch)
        new = flip(sd(toswitch(i),:));
        sd(toswitch(i),:) = new;
    end
end


if size(mu,2) == 1 % if only 1 null dim, squeeze() tranposes for some reason
    mu = mu';
    sd = sd';
end
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.null.cdEarly_mode = cd;

% potent
tempdat = rez.potent_psth;
mu = squeeze(mean(tempdat(times.early,:,:),1));
if switching
    toswitch = find(mu(:,1) < mu(:,2));
    for i = 1:numel(toswitch)
        new = flip(mu(toswitch(i),:));
        mu(toswitch(i),:) = new;
    end
end

sd = squeeze(std(tempdat(times.early,:,:),[],1));
if switching
    for i = 1:numel(toswitch)
        new = flip(sd(toswitch(i),:));
        sd(toswitch(i),:) = new;
    end
end


if size(mu,2) == 1 % if only 1 null dim, squeeze() tranposes for some reason
    mu = mu';
    sd = sd';
end
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.potent.cdEarly_mode = cd;

%% cd late

% null
tempdat = rez.null_psth;
mu = squeeze(mean(tempdat(times.late,:,:),1));
if switching
    toswitch = find(mu(:,1) < mu(:,2));
    for i = 1:numel(toswitch)
        new = flip(mu(toswitch(i),:));
        mu(toswitch(i),:) = new;
    end
end

sd = squeeze(std(tempdat(times.late,:,:),[],1));
if switching
    for i = 1:numel(toswitch)
        new = flip(sd(toswitch(i),:));
        sd(toswitch(i),:) = new;
    end
end


if size(mu,2) == 1 % if only 1 null dim, squeeze() tranposes for some reason
    mu = mu';
    sd = sd';
end
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.null.cdLate_mode = cd;

% potent
tempdat = rez.potent_psth;
mu = squeeze(mean(tempdat(times.late,:,:),1));
if switching
    toswitch = find(mu(:,1) < mu(:,2));
    for i = 1:numel(toswitch)
        new = flip(mu(toswitch(i),:));
        mu(toswitch(i),:) = new;
    end
end

sd = squeeze(std(tempdat(times.late,:,:),[],1));
if switching
    for i = 1:numel(toswitch)
        new = flip(sd(toswitch(i),:));
        sd(toswitch(i),:) = new;
    end
end


if size(mu,2) == 1 % if only 1 null dim, squeeze() tranposes for some reason
    mu = mu';
    sd = sd';
end
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.potent.cdLate_mode = cd;

%% cd go

% null
tempdat = rez.null_psth;
mu = squeeze(mean(tempdat(times.go,:,:),1));
if switching
    toswitch = find(mu(:,1) < mu(:,2));
    for i = 1:numel(toswitch)
        new = flip(mu(toswitch(i),:));
        mu(toswitch(i),:) = new;
    end
end

sd = squeeze(std(tempdat(times.go,:,:),[],1));
if switching
    for i = 1:numel(toswitch)
        new = flip(sd(toswitch(i),:));
        sd(toswitch(i),:) = new;
    end
end


if size(mu,2) == 1 % if only 1 null dim, squeeze() tranposes for some reason
    mu = mu';
    sd = sd';
end
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.null.cdGo_mode = cd;

% potent
tempdat = rez.potent_psth;
mu = squeeze(mean(tempdat(times.go,:,:),1));
if switching
    toswitch = find(mu(:,1) < mu(:,2));
    for i = 1:numel(toswitch)
        new = flip(mu(toswitch(i),:));
        mu(toswitch(i),:) = new;
    end
end

sd = squeeze(std(tempdat(times.go,:,:),[],1));
if switching
    for i = 1:numel(toswitch)
        new = flip(sd(toswitch(i),:));
        sd(toswitch(i),:) = new;
    end
end


if size(mu,2) == 1 % if only 1 null dim, squeeze() tranposes for some reason
    mu = mu';
    sd = sd';
end
cd = ((mu(:,1)-mu(:,2)))./ sqrt(sum(sd.^2,2));
cd(isnan(cd)) = 0;
cd = cd./sum(abs(cd)); % (ncells,1)
rez.cd.potent.cdGo_mode = cd;

%% orthogonalize

[fns,~] = patternMatchCellArray(fieldnames(rez.cd.null),{'mode'},'all');
nullmodes = zeros(rez.dPrep,numel(fns));
potentmodes = zeros(rez.dMove,numel(fns));
for i = 1:numel(fns)
    nullmodes(:,i) = rez.cd.null.(fns{i});
    potentmodes(:,i) = rez.cd.potent.(fns{i});
end

orthNullModes = gschmidt(nullmodes);
orthPotentModes = gschmidt(potentmodes);

for i = 1:numel(fns)
    if orthog
        rez.cd.null.(fns{i}) = orthNullModes(:,i);
        rez.cd.potent.(fns{i}) = orthPotentModes(:,i);
    else
        rez.cd.null.(fns{i}) = nullmodes(:,i);
        rez.cd.potent.(fns{i}) = potentmodes(:,i);
    end
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



end

















