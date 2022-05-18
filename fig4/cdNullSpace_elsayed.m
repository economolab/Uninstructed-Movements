function [latents,trials_by_type,modes] = cdNullSpace_elsayed(rez,dat,obj,params)


%% trial types

trials.rhit = 1:numel(rez.rhit);
trials.lhit = numel(rez.rhit)+1:numel(rez.rhit)+numel(rez.lhit);

%% epochs

alignment = mode(obj.bp.ev.(params.alignEvent));

% early sample epoch
e1 = obj.bp.ev.sample(dat.trials) - alignment; 
e2 = obj.bp.ev.sample(dat.trials) + 0.4 - alignment;
early = obj.time>=mode(e1) & obj.time <= mode(e2);

% late delay epoch
e1 = obj.bp.ev.goCue(dat.trials) - 0.41 - alignment; 
e2 = obj.bp.ev.goCue(dat.trials) - 0.01 - alignment;
late = obj.time>=mode(e1) & obj.time <= mode(e2);

% go epoch
e1 = obj.bp.ev.goCue(dat.trials) + 0.01 - alignment; 
e2 = obj.bp.ev.goCue(dat.trials) + 0.41 - alignment;
go = obj.time>=mode(e1) & obj.time <= mode(e2);

%% early sample cd

conditions = fieldnames(trials);
datnull = cell(numel(conditions),1);
datpotent = cell(numel(conditions),1);
meannull_time = cell(numel(conditions),1);
meanpotent_time = cell(numel(conditions),1);
mu_null = nan(size(rez.N_null,2),numel(conditions));
mu_potent = nan(size(rez.N_potent,2),numel(conditions));
sd_null = nan(size(rez.N_null,2),numel(conditions));
sd_potent = nan(size(rez.N_potent,2),numel(conditions));
for i = 1:numel(conditions)
    datnull{i} = rez.N_null(early,:,trials.(conditions{i}));
    datpotent{i} = rez.N_potent(early,:,trials.(conditions{i}));
    
    meannull_time{i} = squeeze(mean(datnull{i},1));
    meanpotent_time{i} = squeeze(mean(datpotent{i},1));
    
    mu_null(:,i) = nanmean(meannull_time{i},2);
    mu_potent(:,i) = nanmean(meanpotent_time{i},2);
   
    sd_null(:,i) = nanstd(meannull_time{i},[],2);
    sd_potent(:,i) = nanstd(meanpotent_time{i},[],2);
    
end


% calculate mode according to definition
modes.null.early = (mu_null(:,1)-mu_null(:,2))./ sqrt(sum(sd_null.^2,2));
modes.null.early(isnan(modes.null.early)) = 0;
modes.null.early = modes.null.early./sum(abs(modes.null.early)); % (ncells,1)

modes.potent.early = (mu_potent(:,1)-mu_potent(:,2)) ./ sqrt(sum(sd_potent.^2,2));
modes.potent.early(isnan(modes.potent.early)) = 0;
modes.potent.early = modes.potent.early./sum(abs(modes.potent.early)); % (ncells,1)


%% late delay cd

conditions = fieldnames(trials);
datnull = cell(numel(conditions),1);
datpotent = cell(numel(conditions),1);
meannull_time = cell(numel(conditions),1);
meanpotent_time = cell(numel(conditions),1);
mu_null = nan(size(rez.N_null,2),numel(conditions));
mu_potent = nan(size(rez.N_potent,2),numel(conditions));
sd_null = nan(size(rez.N_null,2),numel(conditions));
sd_potent = nan(size(rez.N_potent,2),numel(conditions));
for i = 1:numel(conditions)
    datnull{i} = rez.N_null(late,:,trials.(conditions{i}));
    datpotent{i} = rez.N_potent(late,:,trials.(conditions{i}));
    
    meannull_time{i} = squeeze(mean(datnull{i},1));
    meanpotent_time{i} = squeeze(mean(datpotent{i},1));
    
    mu_null(:,i) = nanmean(meannull_time{i},2);
    mu_potent(:,i) = nanmean(meanpotent_time{i},2);
   
    sd_null(:,i) = nanstd(meannull_time{i},[],2);
    sd_potent(:,i) = nanstd(meanpotent_time{i},[],2);
    
end


% calculate mode according to definition
modes.null.late = (mu_null(:,1)-mu_null(:,2))./ sqrt(sum(sd_null.^2,2));
modes.null.late(isnan(modes.null.late)) = 0;
modes.null.late = modes.null.late./sum(abs(modes.null.late)); % (ncells,1)

modes.potent.late = (mu_potent(:,1)-mu_potent(:,2)) ./ sqrt(sum(sd_potent.^2,2));
modes.potent.late(isnan(modes.potent.late)) = 0;
modes.potent.late = modes.potent.late./sum(abs(modes.potent.late)); % (ncells,1)

%% go mode

conditions = fieldnames(trials);
datnull = cell(numel(conditions),1);
datpotent = cell(numel(conditions),1);
meannull_time = cell(numel(conditions),1);
meanpotent_time = cell(numel(conditions),1);
mu_null = nan(size(rez.N_null,2),numel(conditions));
mu_potent = nan(size(rez.N_potent,2),numel(conditions));
sd_null = nan(size(rez.N_null,2),numel(conditions));
sd_potent = nan(size(rez.N_potent,2),numel(conditions));
for i = 1:numel(conditions)
    datnull{i} = rez.N_null(go,:,trials.(conditions{i}));
    datpotent{i} = rez.N_potent(go,:,trials.(conditions{i}));
    
    meannull_time{i} = squeeze(mean(datnull{i},1));
    meanpotent_time{i} = squeeze(mean(datpotent{i},1));
    
    mu_null(:,i) = nanmean(meannull_time{i},2);
    mu_potent(:,i) = nanmean(meanpotent_time{i},2);
   
    sd_null(:,i) = nanstd(meannull_time{i},[],2);
    sd_potent(:,i) = nanstd(meanpotent_time{i},[],2);
    
end

% calcugo mode according to definition
modes.null.go = (mu_null(:,1)-mu_null(:,2))./ sqrt(sum(sd_null.^2,2));
modes.null.go(isnan(modes.null.go)) = 0;
modes.null.go = modes.null.go./sum(abs(modes.null.go)); % (ncells,1)

modes.potent.go = (mu_potent(:,1)-mu_potent(:,2)) ./ sqrt(sum(sd_potent.^2,2));
modes.potent.go(isnan(modes.potent.go)) = 0;
modes.potent.go = modes.potent.go./sum(abs(modes.potent.go)); % (ncells,1)

%% othogonalize

nullModes = [modes.null.early modes.null.late modes.null.go];
nullModes_orth = gschmidt(nullModes);

potentModes = [modes.potent.early modes.potent.late modes.potent.go];
potentModes_orth = gschmidt(potentModes);

%% projections

% plot null/potent projected neural data onto null activity modes
proj_null_null = zeros(size(rez.N_null,1),size(nullModes_orth,2),size(rez.N_null,3)); % projection of rez.N_null onto activity modes defined by null space
% proj_potent_null = zeros(size(rez.N_null,1),size(nullModes_orth,2),size(rez.N_potent,3)); % projection of rez.N_potent onto activity modes defined by null space
for trix = 1:size(rez.N_null,3) % number of trials
    temp = squeeze(rez.N_null(:,:,trix)); % (time,'numClusters')
    proj_null_null(:,:,trix) = temp * nullModes_orth;
%     temp = squeeze(rez.N_potent(:,:,trix));
%     proj_potent_null(:,:,trix) = temp * nullModes_orth;
end

% plot null/potent projected neural data onto potent activity modes
% proj_null_potent = zeros(size(rez.N_null,1),size(potentModes_orth,2),size(rez.N_null,3)); % projection of rez.N_null onto activity modes defined by null space
proj_potent_potent = zeros(size(rez.N_null,1),size(potentModes_orth,2),size(rez.N_potent,3)); % projection of rez.N_potent onto activity modes defined by null space
for trix = 1:size(rez.N_null,3) % number of trials
%     temp = squeeze(rez.N_null(:,:,trix)); % (time,'numClusters')
%     proj_null_potent(:,:,trix) = temp * potentModes_orth;
    temp = squeeze(rez.N_potent(:,:,trix));
    proj_potent_potent(:,:,trix) = temp * potentModes_orth;
end

%% clean up

fns = fieldnames(modes.null);

latents.names = fns;
latents.proj_null_null = proj_null_null;
% latents.proj_potent_null = proj_potent_null;
% latents.proj_null_potent = proj_null_potent;
latents.proj_potent_potent = proj_potent_potent;
latents.time = obj.time;
latents.plotTimeIx = 15:numel(obj.time);

trials_by_type.rhit = trials.rhit;
trials_by_type.lhit = trials.lhit;

end

















