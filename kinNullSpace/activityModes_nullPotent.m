function [latents,trials_by_type,modes] = activityModes_nullPotent(rez,dat,obj,params)


%% trial types

rhit = find(obj.bp.R & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
lhit = find(obj.bp.L & obj.bp.hit & ~obj.bp.autowater & ~obj.bp.early);
mask = ismember(dat.trials,rhit);
trials.rhit = find(mask);
mask = ismember(dat.trials,lhit);
trials.lhit = find(mask);

rmiss = find(obj.bp.R & obj.bp.miss & ~obj.bp.autowater & ~obj.bp.early);
lmiss = find(obj.bp.L & obj.bp.miss & ~obj.bp.autowater & ~obj.bp.early);
mask = ismember(dat.trials,rmiss);
trials.rmiss = find(mask);
mask = ismember(dat.trials,lmiss);
trials.lmiss = find(mask);

hits = find(obj.bp.hit & ~obj.bp.autowater & ~ obj.bp.early);
mask = ismember(dat.trials,hits);
hits = find(mask);

%% epochs

alignment = mode(obj.bp.ev.(params.alignEvent));

% presample epoch
e1 = obj.bp.ev.sample(dat.trials) - 0.5 - alignment; 
e2 = obj.bp.ev.sample(dat.trials) - 0.05 - alignment;
presample = obj.time>=mode(e1) & obj.time <= mode(e2);

% sample epoch
e1 = obj.bp.ev.sample(dat.trials) - alignment;
e2 = obj.bp.ev.sample(dat.trials) - 0.3 - alignment;
sample = obj.time>=mode(e2) & obj.time <= mode(e1);

% delay epoch
e1 = obj.bp.ev.delay(dat.trials) - alignment;
e2 = obj.bp.ev.goCue(dat.trials) - 0.05 - alignment;
delay = obj.time>=mode(e1) & obj.time <= mode(e2);

% action epoch
e1 = obj.bp.ev.goCue(dat.trials) + 0.1 - alignment;
e2 = obj.bp.ev.goCue(dat.trials) + 0.3 - alignment;
action = obj.time>=mode(e1) & obj.time <= mode(e2);

% outcome epoch
e1 = obj.bp.ev.goCue(dat.trials) - alignment;
e2 = obj.bp.ev.goCue(dat.trials) + 1.3 - alignment;
outcome = obj.time>=mode(e1) & obj.time <= mode(e2);

% pre go cue epoch
e1 = obj.bp.ev.goCue(dat.trials) - 0.1 - alignment;
e2 = obj.bp.ev.goCue(dat.trials) - alignment;
prego = obj.time>=mode(e1) & obj.time <= mode(e2);

% post go cue epoch
e1 = obj.bp.ev.goCue(dat.trials) - alignment;
e2 = obj.bp.ev.goCue(dat.trials) + 0.1 - alignment;
postgo = obj.time>=mode(e1) & obj.time <= mode(e2);


%% stimulus mode

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
    datnull{i} = rez.N_null(sample,:,trials.(conditions{i}));
    datpotent{i} = rez.N_potent(sample,:,trials.(conditions{i}));
    
    meannull_time{i} = squeeze(mean(datnull{i},1));
    meanpotent_time{i} = squeeze(mean(datpotent{i},1));
    
    mu_null(:,i) = nanmean(meannull_time{i},2);
    mu_potent(:,i) = nanmean(meanpotent_time{i},2);
   
    sd_null(:,i) = nanstd(meannull_time{i},[],2);
    sd_potent(:,i) = nanstd(meanpotent_time{i},[],2);
    
end


% calculate mode according to definition
modes.null.stimulus = ((mu_null(:,1)-mu_null(:,4)) + (mu_null(:,3)-mu_null(:,2)))./ sqrt(sum(sd_null.^2,2));
modes.null.stimulus(isnan(modes.null.stimulus)) = 0;
modes.null.stimulus = modes.null.stimulus./sum(abs(modes.null.stimulus)); % (ncells,1)

modes.potent.stimulus = ((mu_potent(:,1)-mu_potent(:,4)) + (mu_potent(:,3)-mu_potent(:,2)))./ sqrt(sum(sd_potent.^2,2));
modes.potent.stimulus(isnan(modes.potent.stimulus)) = 0;
modes.potent.stimulus = modes.potent.stimulus./sum(abs(modes.potent.stimulus)); % (ncells,1)

%% choice mode

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
    datnull{i} = rez.N_null(delay,:,trials.(conditions{i}));
    datpotent{i} = rez.N_potent(delay,:,trials.(conditions{i}));
    
    meannull_time{i} = squeeze(mean(datnull{i},1));
    meanpotent_time{i} = squeeze(mean(datpotent{i},1));
    
    mu_null(:,i) = nanmean(meannull_time{i},2);
    mu_potent(:,i) = nanmean(meanpotent_time{i},2);
   
    sd_null(:,i) = nanstd(meannull_time{i},[],2);
    sd_potent(:,i) = nanstd(meanpotent_time{i},[],2);
    
end


% calculate mode according to definition
modes.null.choice = ((mu_null(:,1)-mu_null(:,3)) + (mu_null(:,4)-mu_null(:,2)))./ sqrt(sum(sd_null.^2,2));
modes.null.choice(isnan(modes.null.choice)) = 0;
modes.null.choice = modes.null.choice./sum(abs(modes.null.choice)); % (ncells,1)

modes.potent.choice = ((mu_potent(:,1)-mu_potent(:,3)) + (mu_potent(:,4)-mu_potent(:,2)))./ sqrt(sum(sd_potent.^2,2));
modes.potent.choice(isnan(modes.potent.choice)) = 0;
modes.potent.choice = modes.potent.choice./sum(abs(modes.potent.choice)); % (ncells,1)


%% action mode

conditions = fieldnames(trials);
conditions = conditions(1:2);
datnull = cell(numel(conditions),1);
datpotent = cell(numel(conditions),1);
meannull_time = cell(numel(conditions),1);
meanpotent_time = cell(numel(conditions),1);
mu_null = nan(size(rez.N_null,2),numel(conditions));
mu_potent = nan(size(rez.N_potent,2),numel(conditions));
sd_null = nan(size(rez.N_null,2),numel(conditions));
sd_potent = nan(size(rez.N_potent,2),numel(conditions));
for i = 1:numel(conditions)
    datnull{i} = rez.N_null(action,:,trials.(conditions{i}));
    datpotent{i} = rez.N_potent(action,:,trials.(conditions{i}));
    
    meannull_time{i} = squeeze(mean(datnull{i},1));
    meanpotent_time{i} = squeeze(mean(datpotent{i},1));
    
    mu_null(:,i) = nanmean(meannull_time{i},2);
    mu_potent(:,i) = nanmean(meanpotent_time{i},2);
   
    sd_null(:,i) = nanstd(meannull_time{i},[],2);
    sd_potent(:,i) = nanstd(meanpotent_time{i},[],2);
    
end


% calculate mode according to definition
modes.null.action = (mu_null(:,1)-mu_null(:,2)) ./ sqrt(sum(sd_null.^2,2));
modes.null.action(isnan(modes.null.choice)) = 0;
modes.null.action = modes.null.action./sum(abs(modes.null.action)); % (ncells,1)

modes.potent.action = (mu_potent(:,1)-mu_potent(:,2))./ sqrt(sum(sd_potent.^2,2));
modes.potent.action(isnan(modes.potent.action)) = 0;
modes.potent.action = modes.potent.action./sum(abs(modes.potent.action)); % (ncells,1)


%% outcome mode

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
    datnull{i} = rez.N_null(outcome,:,trials.(conditions{i}));
    datpotent{i} = rez.N_potent(outcome,:,trials.(conditions{i}));
    
    meannull_time{i} = squeeze(mean(datnull{i},1));
    meanpotent_time{i} = squeeze(mean(datpotent{i},1));
    
    mu_null(:,i) = nanmean(meannull_time{i},2);
    mu_potent(:,i) = nanmean(meanpotent_time{i},2);
   
    sd_null(:,i) = nanstd(meannull_time{i},[],2);
    sd_potent(:,i) = nanstd(meanpotent_time{i},[],2);
    
end


% calculate mode according to definition
modes.null.outcome = ((mu_null(:,1)-mu_null(:,3)) + (mu_null(:,2)-mu_null(:,4)))./ sqrt(sum(sd_null.^2,2));
modes.null.outcome(isnan(modes.null.outcome)) = 0;
modes.null.outcome = modes.null.outcome./sum(abs(modes.null.outcome)); % (ncells,1)

modes.potent.outcome = ((mu_potent(:,1)-mu_potent(:,3)) + (mu_potent(:,2)-mu_potent(:,4)))./ sqrt(sum(sd_potent.^2,2));
modes.potent.outcome(isnan(modes.potent.outcome)) = 0;
modes.potent.outcome = modes.potent.outcome./sum(abs(modes.potent.outcome)); % (ncells,1)

%% ramping mode

% presample epoch
datnull = rez.N_null(presample,:,hits);
datpotent = rez.N_potent(presample,:,hits);

meannull_time = squeeze(mean(datnull,1));
meanpotent_time = squeeze(mean(datpotent,1));

mu_null_presample = nanmean(meannull_time,2);
mu_potent_presample = nanmean(meanpotent_time,2);

sd_null_presample = nanstd(meannull_time,[],2);
sd_potent_presample = nanstd(meanpotent_time,[],2);

% delay epoch
datnull = rez.N_null(delay,:,hits);
datpotent = rez.N_potent(delay,:,hits);

meannull_time = squeeze(mean(datnull,1));
meanpotent_time = squeeze(mean(datpotent,1));

mu_null_delay = nanmean(meannull_time,2);
mu_potent_delay = nanmean(meanpotent_time,2);

sd_null_delay = nanstd(meannull_time,[],2);
sd_potent_delay = nanstd(meanpotent_time,[],2);


sd_null =  [sd_null_delay sd_null_presample]; 
sd_potent =  [sd_potent_delay sd_potent_presample]; 

% calculate mode according to definition
modes.null.ramping = (mu_null_presample - mu_null_delay) ./ sqrt(sum(sd_null.^2,2));
modes.null.ramping(isnan(modes.null.ramping)) = 0;
modes.null.ramping = modes.null.ramping./sum(abs(modes.null.ramping)); % (ncells,1)

modes.potent.ramping = (mu_potent_presample - mu_potent_delay) ./ sqrt(sum(sd_potent.^2,2));
modes.potent.ramping(isnan(modes.potent.ramping)) = 0;
modes.potent.ramping = modes.potent.ramping./sum(abs(modes.potent.ramping)); % (ncells,1)


%% go mode

% prego epoch
datnull = rez.N_null(prego,:,hits);
datpotent = rez.N_potent(prego,:,hits);

meannull_time = squeeze(mean(datnull,1));
meanpotent_time = squeeze(mean(datpotent,1));

mu_null_prego = nanmean(meannull_time,2);
mu_potent_prego = nanmean(meanpotent_time,2);

sd_null_prego = nanstd(meannull_time,[],2);
sd_potent_prego = nanstd(meanpotent_time,[],2);

% postgo epoch
datnull = rez.N_null(postgo,:,hits);
datpotent = rez.N_potent(postgo,:,hits);

meannull_time = squeeze(mean(datnull,1));
meanpotent_time = squeeze(mean(datpotent,1));

mu_null_postgo = nanmean(meannull_time,2);
mu_potent_postgo = nanmean(meanpotent_time,2);

sd_null_postgo = nanstd(meannull_time,[],2);
sd_potent_postgo = nanstd(meanpotent_time,[],2);


sd_null =  [sd_null_postgo sd_null_prego]; 
sd_potent =  [sd_potent_postgo sd_potent_prego]; 

% calculate mode according to definition
modes.null.go = (mu_null_postgo - mu_null_prego) ./ sqrt(sum(sd_null.^2,2));
modes.null.go(isnan(modes.null.go)) = 0;
modes.null.go = modes.null.go./sum(abs(modes.null.go)); % (ncells,1)

modes.potent.go = (mu_potent_postgo - mu_potent_prego) ./ sqrt(sum(sd_potent.^2,2));
modes.potent.go(isnan(modes.potent.go)) = 0;
modes.potent.go = modes.potent.go./sum(abs(modes.potent.go)); % (ncells,1)

%% othogonalize

nullModes = [modes.null.stimulus modes.null.choice modes.null.action modes.null.outcome modes.null.ramping modes.null.go];
nullModes_orth = gschmidt(nullModes);

potentModes = [modes.potent.stimulus modes.potent.choice modes.potent.action modes.potent.outcome modes.potent.ramping modes.potent.go];
potentModes_orth = gschmidt(potentModes);

%% projections

% plot null/potent projected neural data onto null activity modes
proj_null_null = zeros(size(rez.N_null,1),size(nullModes_orth,2),size(rez.N_null,3)); % projection of rez.N_null onto activity modes defined by null space
proj_potent_null = zeros(size(rez.N_null,1),size(nullModes_orth,2),size(rez.N_potent,3)); % projection of rez.N_potent onto activity modes defined by null space
for trix = 1:size(rez.N_null,3) % number of trials
    temp = squeeze(rez.N_null(:,:,trix)); % (time,'numClusters')
    proj_null_null(:,:,trix) = temp * nullModes_orth;
    temp = squeeze(rez.N_potent(:,:,trix));
    proj_potent_null(:,:,trix) = temp * nullModes_orth;
end

% plot null/potent projected neural data onto potent activity modes
proj_null_potent = zeros(size(rez.N_null,1),size(potentModes_orth,2),size(rez.N_null,3)); % projection of rez.N_null onto activity modes defined by null space
proj_potent_potent = zeros(size(rez.N_null,1),size(potentModes_orth,2),size(rez.N_potent,3)); % projection of rez.N_potent onto activity modes defined by null space
for trix = 1:size(rez.N_null,3) % number of trials
    temp = squeeze(rez.N_null(:,:,trix)); % (time,'numClusters')
    proj_null_potent(:,:,trix) = temp * potentModes_orth;
    temp = squeeze(rez.N_potent(:,:,trix));
    proj_potent_potent(:,:,trix) = temp * potentModes_orth;
end

%% clean up

fns = fieldnames(modes.null);

latents.names = fns;
latents.proj_null_null = proj_null_null;
latents.proj_potent_null = proj_potent_null;
latents.proj_null_potent = proj_null_potent;
latents.proj_potent_potent = proj_potent_potent;
latents.time = obj.time;
latents.plotTimeIx = 15:numel(obj.time);

trials_by_type.rhit = trials.rhit;
trials_by_type.lhit = trials.lhit;
trials_by_type.rmiss = trials.rmiss;
trials_by_type.lmiss = trials.lmiss;
trials_by_type.hit = hits;

end

















